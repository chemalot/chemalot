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
import java.io.IOException;
import java.io.InputStreamReader;

import openeye.oechem.*;

import org.apache.commons.cli.*;


/**

 * @author Ben Sellers / Alberto Gobbi / 2014
 * Copyright 2014 Genentech
 *
 * Class to support generation of torsion conformers from command line
 */
public class SDFTorsionScanner
{
   private static final String MY_NAME = "SDFTorsionScanner";
   public  static final int FRAGNumTag = oechem.OEGetTag("fragNum");

   private static final String OPT_INFILE            = "in";
   private static final String OPT_OUTFILE           = "out";
   private static final String OPT_BONDFILE          = "bondFile";
   private static final String OPT_STARTTorsion      = "startTorsion";
   private static final String OPT_TORSIONIncrement  = "torsionIncrement";
   private static final String OPT_NSTEPS            = "nSteps";
   private static final String OPT_COREFILE          = "core";
   private static final String OPT_TORSIONS_ATOM_TAG = "torsionAtomsTag";
   private static final String OPT_MINIMIZE          = "minimize";
   private static final String OPT_MAXCONFS_PER_STEP = "maxConfsPerStep";
   private static final String OPT_CONSTRIANT        = "constraintStrength";

   private static final String TORMOLNUM_TAG         = "torMolNum";

   private final TorsionScanner torsionScanner;



   /**
    * Constructor
    * @param outFile the file to write to
    * @param bondFile sdf file with four atoms whose positions define the bond in the input file
    * @param coreFilename core filename to write out
    * @param torsionAtomsTag the tag name to which to write out atom indices for the torsion
    * @param nSteps number of angular steps to generate around torsion
    * @param nStartTor the starting torsion angle in degrees
    * @param nTorIncr the degree increment between generated torsion angles
    * @param nMaxConfsPerStep the max number of conformers to generate using rotatable bonds other than the fixed torsion
    * @param doMinimize minimize at each step.  if the nMaxConfsPerStep is > 1, write out only the min E conformer per step
    * @param constraintStrength
    */
   public SDFTorsionScanner(String bondFilename,
                            String coreFilename,
                            String torsionAtomsTag,
                            int nSteps,
                            int nStartTor,
                            int nTorIncr,
                            int nMaxConfsPerStep,
                            boolean doMinimize,
                            String constraintStrength)
   {
      this.torsionScanner   =
        new TorsionScanner( bondFilename,
                            coreFilename,
                            torsionAtomsTag,
                            nSteps,
                            nStartTor,
                            nTorIncr,
                            nMaxConfsPerStep,
                            doMinimize,
                            constraintStrength);
   }


   private void run(String inFile, String outFile, String coreFile) throws IOException
   {

      oemolothread outputOEThread  = new oemolothread(outFile);
      oemolithread inputOEThread   = new oemolithread(inFile);

      OEGraphMol mol = new OEGraphMol();

      Integer nMolIndex = 0;

      try
      {
         // Read each molecule in input
         while (oechem.OEReadMolecule(inputOEThread, mol))
         {
            int     nTotalConformers = 0;

            setMolIndex(mol, nMolIndex);

            //
            // Generate conformers
            OEMCMolBase torsionConformers = torsionScanner.run(mol);
            nTotalConformers = torsionConformers.NumConfs();

            //
            // Write out conformers
            oechem.OEWriteMolecule(outputOEThread, torsionConformers);

            // Cleanup
            torsionConformers.delete();

            System.err.println(torsionScanner.getSteps() + " conformers generated at " + torsionScanner.getTorIncr() + " degree increments.");
            if (torsionScanner.getMaxConfsPerStep() > 1 && !torsionScanner.doMinimize())
            {
               System.err.println("Expanded to " + Integer.toString(nTotalConformers) + " total conformers.");
            }
            nMolIndex++;
         }

         System.err.println(nMolIndex.toString() + " molecule(s) read from input file.");

         if(coreFile != null)
            torsionScanner.computeCore(coreFile);

         mol.delete();
      }
      catch(IllegalArgumentException iaEx)
      {
         System.err.println("Problem computing core file: " + iaEx.getMessage());
      }
      catch(Exception e)
      {
         System.err.println(e.getMessage());
      }
      finally
      {
         torsionScanner.close();
         inputOEThread.close();
         inputOEThread.delete();
         outputOEThread.close();
         outputOEThread.delete();
      }
   }


   /**
    * Set the title and index on the molecule
    * @param mol
    * @param nMolIndex
    */
   private static void setMolIndex(OEGraphMol mol, Integer nMolIndex)
   {
      mol.SetTitle("Molecule_" + nMolIndex.toString());

      // Set tor mol num sd tag
      oechem.OESetSDData(mol, TORMOLNUM_TAG, nMolIndex.toString());
   }


   /**
    * Help
    * @param options
    */
   private static void exitWithHelp(Options options)
   {  HelpFormatter formatter = new HelpFormatter();
      String head = "Generates an output file containing multiple torsion conformations.";
      formatter.printHelp( MY_NAME, head, options, "", true );
      System.exit(1);
   }


   /**
    * @param args
    */
   public static void main( String...args ) throws IOException
   {
      // create command line Options object
      Options options = new Options();
      Option opt = new Option( OPT_INFILE, true,
               "input file oe-supported Use .sdf|.smi to specify the file type." );
      opt.setRequired( true );
      options.addOption( opt );

      opt = new Option(OPT_OUTFILE, true,
               "Output filename.");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option(OPT_STARTTorsion, true,
               "The torsion in your inMol will be rotated by this value for the first job");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option(OPT_TORSIONIncrement, true,
               "Incremnt each subsequent conformation by this step size");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option(OPT_NSTEPS, true, "Number of conformations to create");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option(
               OPT_BONDFILE,
               true,
               "The file containing the bond atoms that define the torsion.");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option(OPT_MINIMIZE, false, "Minimize conformer at each step using MMFFs.  If maxConfsPerStep is > 1, "
                                            + "all confs will be mimimized and the lowest E will be output.");
      opt.setRequired(false);
      options.addOption(opt);


      opt = new Option(OPT_CONSTRIANT, true, "One of strong (90), medium (45), weak(20), none or a floating point number"
                                            +" to specify the tethered constraint strength for -minimize (def=strong)");
      opt.setRequired(false);
      options.addOption(opt);


      opt = new Option(OPT_MAXCONFS_PER_STEP, true, "While holding the torsion fixed, maximum number of conformations of free atoms to generat.  default=1");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option(OPT_COREFILE, true, "Outputfile to store guessed core.");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option(OPT_TORSIONS_ATOM_TAG, true, "Name of sdf tag which will contain the indices of atoms that define the torsion.");
      opt.setRequired(true);
      options.addOption(opt);

      CommandLineParser parser = new PosixParser();
      CommandLine cmd = null;
      try
      {
         cmd = parser.parse( options, args );
      }
      catch( Exception e )
      {
         System.err.println( e.getMessage() );
         exitWithHelp( options );
      }
      args = cmd.getArgs();

      if (args.length != 0)
      {  System.err.println("Unknown arguments" + args);
         exitWithHelp(options);
      }

      if (cmd.hasOption("d"))
      {  System.err.println("Start debugger and press return:");
         new BufferedReader(new InputStreamReader(System.in)).readLine();
      }

      String inFile    = cmd.getOptionValue( OPT_INFILE );
      String outFile   = cmd.getOptionValue( OPT_OUTFILE );

      String bondFile  = cmd.getOptionValue( OPT_BONDFILE );
      String coreFilename = cmd.getOptionValue( OPT_COREFILE );
      String torsionAtomsTag = cmd.getOptionValue( OPT_TORSIONS_ATOM_TAG );

      int nSteps    = Integer.parseInt(cmd.getOptionValue(OPT_NSTEPS));
      int nStartTor = Integer.parseInt(cmd.getOptionValue(OPT_STARTTorsion));
      int nTorIncr  = Integer.parseInt(cmd.getOptionValue(OPT_TORSIONIncrement));

      int nMaxConfsPerStep = 1;
      if (cmd.hasOption(OPT_MAXCONFS_PER_STEP))
      {
         nMaxConfsPerStep = Integer.parseInt(cmd.getOptionValue(OPT_MAXCONFS_PER_STEP));
      }

      String constraintStrength = cmd.getOptionValue(OPT_CONSTRIANT);

      boolean doMinimize=cmd.hasOption(OPT_MINIMIZE);

      SDFTorsionScanner torGenerator = new SDFTorsionScanner
                                                       (bondFile, coreFilename, torsionAtomsTag,
                                                        nSteps, nStartTor, nTorIncr, nMaxConfsPerStep, doMinimize, constraintStrength);

      torGenerator.run(inFile, outFile, coreFilename);
   }

}



