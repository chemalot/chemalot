/*
   Copyright 2008-2015 Genentech Inc.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

*/
package com.genentech.chemistry.openEye.apps;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

import openeye.oechem.*;

import org.apache.commons.cli.*;

import com.aestel.utility.DataFormat;

/**
 *
 * @author Alberto Gobbi/ 2012 Copyright 2012 Genentech
 */
public class SdfRMSDSphereExclusion
{  private static final String MY_NAME = "SdfRMSDSphereExclusion";
   private static OEIsChiralAtom CAF = new OEIsChiralAtom();

   private static final String OPT_INFILE  = "in";
   private static final String OPT_OUTFILE = "out";
   private static final String OPT_RADIUS  = "radius";
   private static final String OPT_REFFILE = "refFile";
   private static final String OPT_PRINT_All="printAll";
   private static final String OPT_MIRROR  = "mirror";
   private static final String OPT_GROUPBY = "groupBy";
   private static final String OPT_DONotOpt= "doNotOptimize";
   private static final String OPT_MAXDEVIATION = "useMaxDeviation";

   private final oemolothread outputOEThread;
   private final double radius;
   private final boolean printAll;
   private final List<OEMolBase> selected = new ArrayList<OEMolBase>();

   private final boolean doMirror;
   private final String  groupByTag;
   private final boolean doOptimize;
   private final boolean useMaxDeviation;

   private SdfRMSDSphereExclusion(String refFile, String outFile, double radius,
               boolean printAll, boolean doMirror, boolean doOptimize, boolean useMaxDeviation, String groupBy)
   {
      this.radius = radius;
      this.printAll = printAll;
      this.doMirror = doMirror;
      this.groupByTag = groupBy;
      this.doOptimize = doOptimize;
      this.useMaxDeviation = useMaxDeviation;



      if( refFile != null && groupBy != null)
         throw new Error("groupBy and refFile may not be used together!");

      if (refFile != null)
      {  OEMolBase mol = new OEGraphMol();
         oemolithread ifs = new oemolithread(refFile);
         while (oechem.OEReadMolecule(ifs, mol))
         {  selected.add(new OEGraphMol(mol));
         }
         ifs.close();
         ifs.delete();
         mol.delete();
      }

      outputOEThread = new oemolothread(outFile);
   }


   private void run( String inFile )
   {  oemolithread ifs = new oemolithread(inFile);
      long start = System.currentTimeMillis();
      int iCounter = 0; //Structures in the SD file.
      int oCounter = 0;
      String lastGroupByValue = "Gliberich not found in file";

      OEMolBase mol = new OEGraphMol();

      inLoop:while( oechem.OEReadMolecule( ifs, mol ) )
      {  iCounter++;
         //Output "." to show that the program is running.
         if( iCounter % 100 == 0 )
            System.err.print(".");
         if( iCounter % 4000 == 0 )
            System.err.printf( " %d %dsec\n", iCounter, (System.currentTimeMillis()-start)/1000);

         double minRMSD = Double.MAX_VALUE;

         if( groupByTag != null )
         {  // check if a new group starts here and clear sleected if that is the case
            String groupByVal = oechem.OEGetSDData(mol, groupByTag);
            if( ! groupByVal.equals(lastGroupByValue) )
               clearSelected();
            lastGroupByValue = groupByVal;
         }

         OEMolBase mirMol = null;
         if( doMirror && selected.size() == 0 && ! isChiral(mol))
         {  mirMol=  new OEGraphMol(mol);
            createMirror(mirMol);
         }



         // loop over already selected centroids and  check if mol is inside sphere
         for(int i=selected.size()-1; i>=0; i--)
         {
            OEMolBase currentMol = selected.get(i);

            double rmsd;
            if (! useMaxDeviation)
               rmsd = oechem.OERMSD(mol, currentMol, true, true, doOptimize);
            else
               rmsd = computeMaxDeviation(mol, currentMol);


            if( rmsd < 0D ) System.err.println("OERMSD returned -1 are you comparing two different structures?");

            if( mirMol != null )
            {
               double mirRmsd;
               if (!useMaxDeviation)
                  mirRmsd = oechem.OERMSD(mirMol, currentMol, true, true, doOptimize);
               else
                  mirRmsd = computeMaxDeviation(mirMol, currentMol);

               if( mirRmsd < rmsd ) rmsd = mirRmsd;
            }

            if( rmsd < radius ) // mol is inside sphere
            {  if( printAll )
               {  oechem.OESetSDData(mol, "sphereIdx", Integer.toString(i) );
                  oechem.OESetSDData(mol, "centroidRMSD", DataFormat.formatNumber(rmsd, "si3") );
                  oechem.OEWriteMolecule(outputOEThread, mol);
               }
               continue inLoop;
            }

            // mol is not inside sphere make it a new cnetorid
            if( rmsd < minRMSD ) minRMSD = rmsd;
         }

         oechem.OESetSDData(mol, "sphereIdx", Integer.toString(selected.size()) );
         oechem.OESetSDData(mol, "includeIdx", Integer.toString(selected.size()) );
         oechem.OESetSDData(mol, "centroidRMSD", "0" );

         selected.add(new OEGraphMol(mol));
         oechem.OEWriteMolecule(outputOEThread, mol);

         if( mirMol != null) mirMol.delete();

         oCounter ++;
      }

      mol.delete();
      ifs.close();
      ifs.delete();
      inFile = inFile.replaceAll( ".*" + Pattern.quote(File.separator), "" );
      System.err.printf( "%s: Read %d structures from %s. Written %d centroids in %d sec\n",
               MY_NAME, iCounter, inFile, oCounter, (System.currentTimeMillis()-start)/1000 );
   }


   public static boolean isChiral( OEMolBase mol )
   {  OEAtomBaseIter it = mol.GetAtoms(CAF);
      while(it.hasNext())
      {  OEAtomBase at = it.next();
         if( at.GetAtomicNum() == OEElemNo.N && at.GetDegree() != 4 )
            continue; // ignore stereo on N

         it.delete();
         return true;
      }
      it.delete();
      return false;
   }


   public double computeMaxDeviation(OEMolBase mol1, OEMolBase mol2)
   {  // make deepcopies of mol1 and mol2
      OEMolBase m1 = new OEGraphMol(mol1);
      OEMolBase m2 = new OEGraphMol(mol2);

      double rmat[] = new double[9];
      double trans[] = new double[3];

      @SuppressWarnings("unused")
      double rmsd = oechem.OERMSD(m1, m2, true, true, doOptimize, rmat, trans);
      oechem.OERotate(m2, rmat);
      oechem.OETranslate(m2, trans);

      hideCHydrogens(m1);  // suppress C-H hydrogens
      hideCHydrogens(m2);
      m1.Sweep();  // collect garbage otherwise it would be buggy.
      m2.Sweep();

      double coords1[] = new double[m1.NumAtoms() * 3];
      double coords2[] = new double[m2.NumAtoms() * 3];
      m1.GetCoords(coords1);
      m2.GetCoords(coords2);

      double maxDev = 0.0;
      double matchPairMaxDev, xDev, yDev, zDev;
      double minMatchPairRMSD;
      double matchPairRMSD;
      int Idx1, Idx2;

      OEQMol qMol1 = new OEQMol(m1);
      qMol1.BuildExpressions(OEExprOpts.DefaultAtoms, OEExprOpts.DefaultBonds);
      OESubSearch newSearch = new OESubSearch(qMol1);
      oechem.OEPrepareSearch(m2, newSearch);
      OEMatchBaseIter matchIter = newSearch.Match(m2);

      while (matchIter.hasNext())   // loop over matching searches
      {  OEMatchBase pair = matchIter.next();
         OEMatchPairAtomIter iter = pair.GetAtoms();
         // for every pair:
         matchPairRMSD   = 0.0;
         matchPairMaxDev = 0.0;
         minMatchPairRMSD = Double.MAX_VALUE;
         while (iter.hasNext())    // loop over matching atom-pair
         {  // pattern is original; target is the arg taken from OESubSearch.Match(arg)
            OEMatchPairAtom atomPair = iter.next();
            Idx1 = atomPair.getPattern().GetIdx(); // from mol1
            Idx2 = atomPair.getTarget().GetIdx();  // from mol2
            xDev = Math.abs(coords1[Idx1*3]  -coords2[Idx2*3]);
            yDev = Math.abs(coords1[Idx1*3+1]-coords2[Idx2*3+1]);
            zDev = Math.abs(coords1[Idx1*3+2]-coords2[Idx2*3+2]);
            double atomDev = Math.pow(xDev, 2) + Math.pow(yDev, 2) + Math.pow(zDev, 2);
            matchPairRMSD += atomDev;
            double atomDistance = Math.sqrt(atomDev);
            if( atomDistance > matchPairMaxDev )
            {  matchPairMaxDev = atomDistance;  }
         }
         matchPairRMSD = matchPairRMSD / m1.NumAtoms();  // divided by num of atoms
         matchPairRMSD = Math.sqrt(matchPairRMSD);                // take square root; now we have real RMSD

         if (matchPairRMSD < minMatchPairRMSD)
         {  minMatchPairRMSD = matchPairRMSD;
            maxDev = matchPairMaxDev;
         }
         iter.delete();  // destructor
      }

      // calling destructors
      matchIter.delete();
      newSearch.delete();
      qMol1.delete();
      m1.delete();
      m2.delete();

      return maxDev;
   }


   private static void hideCHydrogens(OEMolBase m1)
   {
      OEAtomBaseIter atIt = m1.GetAtoms();
      while(atIt.hasNext())
      {  OEAtomBase at = atIt.next();
         if( at.GetAtomicNum() == 6 ) oechem.OESuppressHydrogens(at);
      }
      atIt.delete();
   }


   public static OEMolBase createMirror(OEMolBase fitmol)
   {  float coords[] = new float[fitmol.GetMaxAtomIdx() * 3];
      fitmol.GetCoords(coords);
      for( int i=0; i < coords.length; i+=3 )
         coords[i] = -coords[i];
      fitmol.SetCoords(coords);

      return fitmol;
   }


   private void clearSelected()
   {  for(int i=0; i<selected.size(); i++)
         selected.get(i).delete();

      selected.clear();
   }

   private void close()
   {  outputOEThread.close();
      outputOEThread.delete();
      clearSelected();
   }

   private static void exitWithHelp(Options options)
   {  HelpFormatter formatter = new HelpFormatter();
      String head = "Will run a RMSD based sphere exlcusion algorithm with the given radius.\n"
                   +"All conforamtions in the input and reference file must be of the same structure.";
      formatter.printHelp(MY_NAME, head, options, "", true);
      System.exit(1);
   }

   /**
    * @param args
    */
   public static void main(String... args) throws IOException
   { // create command line Options object
      Options options = new Options();
      Option opt = new Option(OPT_INFILE, true,
               "input file oe-supported Use .sdf|.smi to specify the file type.");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option(OPT_OUTFILE, true,
               "output file oe-supported. Use .sdf|.smi to specify the file type.");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option(OPT_REFFILE, true,
               "Reference file of molecules which define pre-existign exclusion spheres.");
      options.addOption(opt);

      opt = new Option(OPT_RADIUS, true, "Radius of exclusion spehre in A2.");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option(OPT_GROUPBY, true, "Group by fieldname, run sphere exclusion for consecutive groups of records with same value for this field.");
      options.addOption(opt);

      opt = new Option(OPT_DONotOpt, false, "If specified the RMSD is computed without trying to optimize the alignment.");
      options.addOption(opt);

      opt = new Option(OPT_PRINT_All, false, "print all molecule, check includeIdx tag");
      options.addOption(opt);

      opt = new Option(OPT_MIRROR, false, "For non-chiral molecules also try mirror image");
      options.addOption(opt);

      opt = new Option(OPT_MAXDEVIATION, false, "Using Max Deviation Method");
      options.addOption(opt);

      CommandLineParser parser = new PosixParser();
      CommandLine cmd = null;
      try
      {
         cmd = parser.parse(options, args);
      } catch (Exception e)
      {
         System.err.println(e.getMessage());
         exitWithHelp(options);
      }
      args = cmd.getArgs();

      if (cmd.hasOption("d"))
      {
         System.err.println("Start debugger and press return:");
         new BufferedReader(new InputStreamReader(System.in)).readLine();
      }

      String inFile = cmd.getOptionValue(OPT_INFILE);
      String outFile = cmd.getOptionValue(OPT_OUTFILE);
      String refFile = cmd.getOptionValue(OPT_REFFILE);
      String groupBy = cmd.getOptionValue(OPT_GROUPBY);
      boolean doOptimize = ! cmd.hasOption(OPT_DONotOpt);
      double radius = Double.parseDouble(cmd.getOptionValue(OPT_RADIUS));
      boolean printAll = cmd.hasOption(OPT_PRINT_All);
      boolean doMirror = cmd.hasOption(OPT_MIRROR);
      boolean useMaxDeviation = cmd.hasOption(OPT_MAXDEVIATION);


      SdfRMSDSphereExclusion sphereExclusion = new SdfRMSDSphereExclusion(
               refFile, outFile, radius, printAll, doMirror, doOptimize, useMaxDeviation, groupBy);

      sphereExclusion.run(inFile);
      sphereExclusion.close();
   }
}
