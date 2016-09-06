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

   private final oemolothread outputOEThread;
   private final double radius;
   private final boolean printAll;
   private final List<OEMolBase> selected = new ArrayList<OEMolBase>();

   private final boolean doMirror;
   private final String  groupByTag;
   private final boolean doOptimize;


   private SdfRMSDSphereExclusion(String refFile, String outFile, double radius,
               boolean printAll, boolean doMirror, boolean doOptimize, String groupBy)
   {
      this.radius = radius;
      this.printAll = printAll;
      this.doMirror = doMirror;
      this.groupByTag = groupBy;
      this.doOptimize = doOptimize;

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

         for(int i=selected.size()-1; i>=0; i--)
         {  double rmsd = oechem.OERMSD(mol, selected.get(i), true, true, doOptimize);
            if( rmsd < 0D ) System.err.println("OERMSD returned -1 are you comparing two different structures?");

            if( mirMol != null )
            {  double mirRmsd = oechem.OERMSD(mirMol, selected.get(i), true, true, doOptimize);
               if( mirRmsd < rmsd ) rmsd = mirRmsd;
            }

            if( rmsd < radius )
            {  if( printAll )
               {  oechem.OESetSDData(mol, "sphereIdx", Integer.toString(i) );
                  oechem.OESetSDData(mol, "centroidRMSD", DataFormat.formatNumber(rmsd, "si3") );
                  oechem.OEWriteMolecule(outputOEThread, mol);
               }
               continue inLoop;
            }

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

      opt = new Option(OPT_PRINT_All,false, "print all molecule, check includeIdx tag");
      options.addOption(opt);

      opt = new Option(OPT_MIRROR,false, "For non-chiral molecules also try mirror image");
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


      SdfRMSDSphereExclusion sphereExclusion = new SdfRMSDSphereExclusion(
               refFile, outFile, radius, printAll, doMirror, doOptimize, groupBy);

      sphereExclusion.run(inFile);
      sphereExclusion.close();
   }
}
