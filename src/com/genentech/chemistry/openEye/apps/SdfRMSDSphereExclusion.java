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

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.PosixParser;

import com.aestel.utility.DataFormat;
import com.genentech.chemistry.openEye.RingSystemExtractor;
import com.genentech.oechem.tools.Atom;

import openeye.oechem.OEAtomBase;
import openeye.oechem.OEAtomBaseIter;
import openeye.oechem.OEElemNo;
import openeye.oechem.OEExprOpts;
import openeye.oechem.OEGraphMol;
import openeye.oechem.OEIsChiralAtom;
import openeye.oechem.OEMatchBase;
import openeye.oechem.OEMatchBaseIter;
import openeye.oechem.OEMatchPairAtom;
import openeye.oechem.OEMatchPairAtomIter;
import openeye.oechem.OEMolBase;
import openeye.oechem.OEQMol;
import openeye.oechem.OESubSearch;
import openeye.oechem.oechem;
import openeye.oechem.oemolithread;
import openeye.oechem.oemolothread;

/**
 *
 * @author Alberto Gobbi/ 2012 Copyright 2012 Genentech
 */
public class SdfRMSDSphereExclusion
{  public enum AtomSubset {ALL, LARGEST_RING, LARGEST_RING_SYSTEM};

   private static final String MY_NAME = "SdfRMSDSphereExclusion";
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
   private static final String OPT_DATAFIELD = "dataSphereFieldName";
   private static final String OPT_DATARADIUS= "dataSphereRadius";
   private static final String OPT_DEBUG = "debug";
   private static final String OPT_REPLACE_OUTPUT_WITH_ALIGNED = "replaceOutputWithAligned";
   private static final String OPT_USE_LARGEST_RING_FOR_RMSD = "useLargestRingForRMSD";
   private static final String OPT_USE_LARGEST_RING_SYSTEM_FOR_RMSD = "useLargestRingSystemForRMSD";
   private static final String OPT_USE_ONE_SIDE_CHAIN_ATOM_FOR_RMSD = "useOneSideChainAtomForRMSD";
   
   private final oemolothread outputOEThread;
   private final double radius;
   private final String dataField;
   private final double dataRadius;
   private final boolean printAll;
   private final List<OEMolBase> centroids = new ArrayList<OEMolBase>();

   private final boolean doMirror;
   private final String  groupByTag;
   private final boolean doOptimize;
   private final boolean useMaxDeviation;
   private final boolean debug;
   private final boolean replaceOutputWithAligned;
   private boolean useOneSideChainAtomOnRing;


   private final AtomSubset atomSubsetForRMSD;



   private SdfRMSDSphereExclusion(String refFile,   String outFile,   double radius,
                                  String dataField, double dataRadius, 
                                  boolean printAll, boolean doMirror, boolean doOptimize,
                                  boolean useMaxDeviation, String groupBy, boolean debug,
                                  boolean replaceOutputWithAligned,
                                  AtomSubset atomSubset,
                                  boolean useOneSideChainAtomOnRing)
   {
      this.radius = radius;
      this.dataField = dataField;
      this.dataRadius = dataRadius;
      this.printAll = printAll;
      this.doMirror = doMirror;
      this.groupByTag = groupBy;
      this.doOptimize = doOptimize;
      this.useMaxDeviation = useMaxDeviation;
      this.debug = debug;
      this.replaceOutputWithAligned = replaceOutputWithAligned;

      this.atomSubsetForRMSD = atomSubset;

      this.useOneSideChainAtomOnRing = useOneSideChainAtomOnRing;

      if( refFile != null && groupBy != null)
         throw new Error("groupBy and refFile may not be used together!");

      if (refFile != null)
      {  OEMolBase mol = new OEGraphMol();
         oemolithread ifs = new oemolithread(refFile);
         while (oechem.OEReadMolecule(ifs, mol))
         {  centroids.add(new OEGraphMol(mol));
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

      OEMolBase currentMol = new OEGraphMol();

      inLoop:while( oechem.OEReadMolecule( ifs, currentMol ) )
      {  iCounter++;
         //Output "." to show that the program is running.
         if( iCounter % 100 == 0 )
            System.err.print(".");
         if( iCounter % 4000 == 0 )
            System.err.printf( " %d %dsec\n", iCounter, (System.currentTimeMillis()-start)/1000);

         if( groupByTag != null )
         {  // check if a new group starts here and clear sleected if that is the case
            String groupByVal = oechem.OEGetSDData(currentMol, groupByTag);
            if( ! groupByVal.equals(lastGroupByValue) )
               clearSelected();
            lastGroupByValue = groupByVal;
         }

         OEMolBase mirMol = null;
         // Create a mirror image of this molecule if requested, this is not the first
         // centroid, and the molecule is achiral.
         if( doMirror && centroids.size() != 0 && ! isChiral(currentMol))
         {  mirMol =  new OEGraphMol(currentMol);
            createMirror(mirMol);
         }

         String currentDataValueString = null;
         double currentDataValue = 0.; 
         if( dataField != null )
         {  currentDataValueString = oechem.OEGetSDData(currentMol, dataField);
            if(currentDataValueString != null ) 
               currentDataValue = Double.parseDouble(currentDataValueString);
         }

         // loop over already selected centroids and  check if mol is inside sphere
         for(int i=centroids.size()-1; i>=0; i--)
         {
            OEMolBase currentCentroid = centroids.get(i);

            // first check dataField Delta because it is faster
            if( dataField != null && currentDataValueString != null)
            {  String centroidDataValueS = oechem.OEGetSDData(currentCentroid, dataField);
               if( centroidDataValueS != null )
               {  double delta = currentDataValue - Double.parseDouble(centroidDataValueS);
                  if( delta > dataRadius )
                     continue;  // delta value is larger than radius, cannot in sphere
               }
            }
            
            double rmsd;
            if (! useMaxDeviation )
               rmsd = computeRMSD(currentCentroid, currentMol);
            else
               rmsd = computeMaxDeviation(currentCentroid, currentMol);

            if( rmsd < 0D ) System.err.println("OERMSD returned -1 are you comparing two different structures?");

            if( mirMol != null && (rmsd >= radius || replaceOutputWithAligned || printAll) )
            {
               double mirRmsd;
               if (! useMaxDeviation )
                  mirRmsd = computeRMSD(currentCentroid, mirMol);
               else
                  mirRmsd = computeMaxDeviation(currentCentroid, mirMol);

               if( mirRmsd < 0D ) System.err.println("OERMSD returned -1 are you comparing two different structures?");

               if( mirRmsd < rmsd )
               {  rmsd = mirRmsd;
                  if( replaceOutputWithAligned )
                  {  OEMolBase dmy = currentMol;
                     currentMol = mirMol;
                     mirMol = dmy;
                  }
               }
            }

            if( rmsd < radius ) // mol is inside sphere
            {  if( printAll )
               {  oechem.OESetSDData(currentMol, "sphereIdx", Integer.toString(i) );
                  oechem.OESetSDData(currentMol, "centroidRMSD", DataFormat.formatNumber(rmsd, "si3") );
                  oechem.OEWriteMolecule(outputOEThread, currentMol);
               }
               continue inLoop;
            }

            // mol is not inside sphere make it a new centroid
         }

         oechem.OESetSDData(currentMol, "sphereIdx", Integer.toString(centroids.size()) );
         oechem.OESetSDData(currentMol, "includeIdx", Integer.toString(centroids.size()) );
         oechem.OESetSDData(currentMol, "centroidRMSD", "0" );

         centroids.add(new OEGraphMol(currentMol));
         oechem.OEWriteMolecule(outputOEThread, currentMol);

         if( mirMol != null) mirMol.delete();

         oCounter ++;
      }

      currentMol.delete();
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


   /**
    * Return a subset of atoms from the input molecule
    * @param mol the input molecule
    * @param subset a designation for the desired subset of atoms
    * @return a molecule containing the desired subset
    */
   private OEMolBase getAtomSubset(OEMolBase mol, AtomSubset subset)
   {
      OEMolBase returnMol = mol;
      RingSystemExtractor ringExtractor = null;

      if (subset == AtomSubset.LARGEST_RING )
      {
         ringExtractor = new RingSystemExtractor(useOneSideChainAtomOnRing);
         ringExtractor.extract(mol);

         // Replace the incoming mol with just the largest rings
         if (ringExtractor.hasLargestRing())
         {  returnMol = ringExtractor.getLargestRingMol(); }
      }
      else if (subset == AtomSubset.LARGEST_RING_SYSTEM )
      {
         ringExtractor = new RingSystemExtractor(useOneSideChainAtomOnRing);
         ringExtractor.extract(mol);

         // Replace the incoming mol with just the largest rings
         if (ringExtractor.hasLargestRingSystem() )
         {  returnMol = ringExtractor.getLargestRingSystemMol(); }
      }
      // Else unchanged if no ring
      return returnMol;
   }


   /**
    * Compute the maximum RMSD between the reference refMol and mol following an optimization.
    *
    * Align by RMSD and return the RMSD.  Coordinates of mol will be modified if
    *    replaceOutputWithAligned is true.
    */
   private double computeRMSD(OEMolBase refMol, OEMolBase mol)
   {
      double rmat [] = null;
      double trans[] = null;

      OEMolBase referenceMol = refMol;
      OEMolBase targetMol    = mol;

      // Get a subset of atoms if desired for the RMSD calculation
      referenceMol = getAtomSubset(refMol, atomSubsetForRMSD);
      targetMol    = getAtomSubset(mol,    atomSubsetForRMSD);

      if(debug)
      {
         System.err.println(oechem.OEMolToSmiles(referenceMol));
         System.err.println(oechem.OEMolToSmiles(targetMol));
      }

      double rmsd = -1.0;
      if ( replaceOutputWithAligned )
      {
         rmat  = new double[9];
         trans = new double[3];
         rmsd = oechem.OERMSD(referenceMol, targetMol, true, true, doOptimize, rmat, trans);

         if( rmsd < 0D )
         {
            // ref and target may not have same atoms
            return rmsd;
         }

         // Apply transformation to full molecule
         oechem.OERotate   (mol, rmat);
         oechem.OETranslate(mol, trans);
      }
      else
      {
         rmsd = oechem.OERMSD(referenceMol, targetMol, true, true, doOptimize);
      }
      return rmsd;
   }

   /**
    * Compute the maximum deviation from mol1 and mol2.
    *
    * First align them by RMSD then loop over all atoms to see which atom deviates most form the other
    */
   public double computeMaxDeviation(OEMolBase refMol, OEMolBase mol)
   {
      OEMolBase referenceMol = refMol;
      OEMolBase targetMol    = mol;

      if ( ! replaceOutputWithAligned )
      {
         // make deep copies of mol1 and mol2.  Leave the input mols untouched.
         referenceMol = new OEGraphMol(refMol);
         targetMol    = new OEGraphMol(mol);
      }

      double rmat[]  = new double[9];
      double trans[] = new double[3];

      // Get a subset of atoms if desired for the RMSD calculation
      OEMolBase refRMSDAtoms  = getAtomSubset(referenceMol, atomSubsetForRMSD);
      OEMolBase targRMSDAtoms = getAtomSubset(targetMol,    atomSubsetForRMSD);

      if(debug)
      {
         System.err.println(oechem.OEMolToSmiles(refRMSDAtoms));
         System.err.println(oechem.OEMolToSmiles(targRMSDAtoms));
      }

      double rmsd = oechem.OERMSD(refRMSDAtoms, targRMSDAtoms, true, true, doOptimize, rmat, trans);

      if( rmsd < 0D )
      {
         // ref and target may not have same atoms
         return rmsd;
      }

      // Perform rotation and translation on atom subset for atom matching
      oechem.OERotate   (targRMSDAtoms, rmat);
      oechem.OETranslate(targRMSDAtoms, trans);

      hideCHydrogens(refRMSDAtoms);  // suppress C-H hydrogens
      hideCHydrogens(targRMSDAtoms);
      refRMSDAtoms.Sweep();  // collect garbage in mol object or oebug will result in mismatches.
      targRMSDAtoms.Sweep();

      double coords1[] = new double[refRMSDAtoms.NumAtoms() * 3];
      double coords2[] = new double[targRMSDAtoms.NumAtoms() * 3];
      refRMSDAtoms.GetCoords(coords1);
      targRMSDAtoms.GetCoords(coords2);

      double maxDev = 0.0;
      double matchPairMaxDev, xDev, yDev, zDev;
      double minMatchPairRMSD = Double.MAX_VALUE;
      double matchPairRMSD;
      int Idx1, Idx2;

      // match mol1 to mol2 so that we can make sure we used the best symmetry equivalent match

      OEQMol qMol1 = new OEQMol(refRMSDAtoms);
      qMol1.BuildExpressions(OEExprOpts.DefaultAtoms, OEExprOpts.DefaultBonds);
      OESubSearch newSearch = new OESubSearch(qMol1);
      oechem.OEPrepareSearch(targRMSDAtoms, newSearch);
      OEMatchBaseIter matchIter = newSearch.Match(targRMSDAtoms);

      while (matchIter.hasNext())   // loop over matching searches
      {  OEMatchBase pair = matchIter.next();
         // for every pair:
         matchPairRMSD   = 0.0;
         matchPairMaxDev = 0.0;
         OEMatchPairAtomIter iter = pair.GetAtoms();
         while (iter.hasNext())    // loop over matching atom-pair
         {  // pattern is original; target is the arg taken from OESubSearch.Match(arg)
            OEMatchPairAtom atomPair = iter.next();
            Idx1 = atomPair.getPattern().GetIdx(); // from mol1
            Idx2 = atomPair.getTarget().GetIdx();  // from mol2

            xDev = coords1[Idx1*3]  -coords2[Idx2*3];
            yDev = coords1[Idx1*3+1]-coords2[Idx2*3+1];
            zDev = coords1[Idx1*3+2]-coords2[Idx2*3+2];

            double atomDev = Math.pow(xDev, 2) + Math.pow(yDev, 2) + Math.pow(zDev, 2);
            matchPairRMSD += atomDev;
            double atomDistance = Math.sqrt(atomDev);
            if( atomDistance > matchPairMaxDev )
            {  matchPairMaxDev = atomDistance;
               if(debug)
                  System.err.printf("New distant atom pair: %s-%s=%f\n",
                           Atom.getAtomName(atomPair.getPattern()), Atom.getAtomName(atomPair.getTarget()),
                           atomDistance);
            }
         }
         matchPairRMSD = matchPairRMSD / refRMSDAtoms.NumAtoms();  // divided by num of atoms
         matchPairRMSD = Math.sqrt(matchPairRMSD);                // take square root; now we have real RMSD

         if (matchPairRMSD < minMatchPairRMSD)
         {  if(debug)
               System.err.printf("MaxDeviation oldRMSD=%g newRMSD=%f maxDeviation=%f\n",
                        minMatchPairRMSD, matchPairRMSD, matchPairMaxDev);

            minMatchPairRMSD = matchPairRMSD;
            maxDev = matchPairMaxDev;
         }
         iter.delete();  // destructor
      }

      // calling destructors
      matchIter.delete();
      newSearch.delete();
      qMol1.delete();

      if ( ! replaceOutputWithAligned )
      {
         referenceMol.delete();
         targetMol.delete();
      }

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
   {  for(int i=0; i<centroids.size(); i++)
         centroids.get(i).delete();

      centroids.clear();
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

      opt = new Option(OPT_DATAFIELD, true, "If given: a record is outside of sphere if either the rmsd is > radius or delta-value is > dataRadius"); 
      options.addOption(opt);

      opt = new Option(OPT_DATARADIUS, true, "cf. " + OPT_DATAFIELD + " radius to use on delta value");
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

      opt = new Option(OPT_DEBUG, false, "Additional logging");
      options.addOption(opt);

      opt = new Option(OPT_REPLACE_OUTPUT_WITH_ALIGNED, false, "If specified, the output molecule coordinates will be replaced  with new coordinates aligned to the centroids.");
      options.addOption(opt);

      opt = new Option(OPT_USE_LARGEST_RING_FOR_RMSD, false, "If specified, all RMSD calculations will be performed on the largest ring.  If no rings, then use the whole molecule.");
      options.addOption(opt);

      opt = new Option(OPT_USE_LARGEST_RING_SYSTEM_FOR_RMSD, false, "If specified, all RMSD calculations will be performed on the largest ring system.  If no rings, then use the whole molecule.");
      options.addOption(opt);

      opt = new Option(OPT_USE_ONE_SIDE_CHAIN_ATOM_FOR_RMSD, false, "Use with one of the ring RMSD options.  Adds one side chain atom to the RMSD calculation.");
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

      if (cmd.hasOption(OPT_USE_LARGEST_RING_FOR_RMSD) && cmd.hasOption(OPT_USE_LARGEST_RING_SYSTEM_FOR_RMSD))
      {
         System.err.println( "Please use either largest ring or largest ring system for the RMSD calculation." );
         exitWithHelp( options );
      }
      if( (cmd.hasOption(OPT_DATAFIELD) && ! cmd.hasOption(OPT_DATARADIUS)) 
            || (! cmd.hasOption(OPT_DATAFIELD) && cmd.hasOption(OPT_DATARADIUS)))
      {  System.err.printf( "%s and %s must be specified together\n", OPT_DATAFIELD, OPT_DATARADIUS );
         exitWithHelp( options );
      }

      double dataRadius = 0.;
      
      String inFile = cmd.getOptionValue(OPT_INFILE);
      String outFile = cmd.getOptionValue(OPT_OUTFILE);
      String refFile = cmd.getOptionValue(OPT_REFFILE);
      String groupBy = cmd.getOptionValue(OPT_GROUPBY);
      boolean doOptimize = ! cmd.hasOption(OPT_DONotOpt);
      double radius = Double.parseDouble(cmd.getOptionValue(OPT_RADIUS));
      String dataField = cmd.getOptionValue(OPT_DATAFIELD);
      if( dataField != null) dataRadius= Double.parseDouble(cmd.getOptionValue(OPT_DATARADIUS));
      boolean printAll = cmd.hasOption(OPT_PRINT_All);
      boolean doMirror = cmd.hasOption(OPT_MIRROR);
      boolean useMaxDeviation = cmd.hasOption(OPT_MAXDEVIATION);
      boolean debug = cmd.hasOption(OPT_DEBUG);
      boolean replaceOutputWithAligned = cmd.hasOption(OPT_REPLACE_OUTPUT_WITH_ALIGNED);
      boolean useLargestRingForRMSD = cmd.hasOption(OPT_USE_LARGEST_RING_FOR_RMSD);
      boolean useLargestRingSystemForRMSD = cmd.hasOption(OPT_USE_LARGEST_RING_SYSTEM_FOR_RMSD);
      boolean useOneSideChainAtomForRMSD = cmd.hasOption(OPT_USE_ONE_SIDE_CHAIN_ATOM_FOR_RMSD);

      AtomSubset atomSubset = AtomSubset.ALL;
      if (useLargestRingForRMSD)
      {  atomSubset = AtomSubset.LARGEST_RING; }
      else if (useLargestRingSystemForRMSD)
      {  atomSubset = AtomSubset.LARGEST_RING_SYSTEM; }

      SdfRMSDSphereExclusion sphereExclusion = new SdfRMSDSphereExclusion(
               refFile, outFile, radius, dataField, dataRadius, printAll, doMirror,
               doOptimize, useMaxDeviation, groupBy, debug, replaceOutputWithAligned,
               atomSubset, useOneSideChainAtomForRMSD);

      sphereExclusion.run(inFile);
      sphereExclusion.close();
   }
}
