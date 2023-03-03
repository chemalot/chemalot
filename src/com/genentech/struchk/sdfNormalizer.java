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
package com.genentech.struchk;

import java.io.File;
import java.net.URL;
import java.util.List;

import openeye.oechem.*;

import org.apache.commons.cli.*;

import com.aestel.utility.Message;
import com.genentech.struchk.oeStruchk.OEStruchk;
import com.genentech.struchk.oeStruchk.StruChkHelper.CHECKConfig;
import com.genentech.struchk.oeStruchk.StructFlagAnalysisInterface;

/**
 * @author albertgo 2008
 *
 * Good tests are: c1[nH]ncn1 and c1[nH]cnn1
 */
public class sdfNormalizer {

   enum OUTMolFormat
   {  ORIGINAL, NORMALIZED, TAUTOMERIC, STEREOPARENT }

   private static final String NON_TETRAHEDRAL_CHIRAL_TAG = "NonTetrahedralChiral";

   public static void main(String [] args) {
      long start = System.currentTimeMillis();
      int nMessages = 0;
      int nErrors = 0;
      int nStruct = 0;

      // create command line Options object
      Options options = new Options();
      Option opt = new Option("in",true, "input file [.ism,.sdf,...]");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("out",true, "output file");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("mol",true, "molFile used for output: ORIGINAL(def)|STEREOPARENT|NORMALIZED|TAUTOMERIC");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("shortMessage",false, "Limit message to first 80 characters to conform with sdf file specs.");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("f",true, "configuration file. (if none the default will be read from the jar file).");
      opt.setRequired(false);
      options.addOption(opt);


      CommandLineParser parser = new PosixParser();
      CommandLine cmd;
      try
      {  cmd = parser.parse(options, args);
      } catch (ParseException e)
      {  exitWithHelp(options, e.getMessage());
         throw new Error(e);  // avoid compiler errors
      }
      args = cmd.getArgs();

      if(args.length != 0) {
         System.err.print("Unknown options: "+ args + "\n\n");
         HelpFormatter formatter = new HelpFormatter();
         formatter.printHelp( "sdfNormalizer", options );
         System.exit(1);
      }

      String molOpt = cmd.getOptionValue("mol");
      OUTMolFormat outMol = OUTMolFormat.ORIGINAL;
      if( molOpt == null || "original".equalsIgnoreCase(molOpt) )
         outMol = OUTMolFormat.ORIGINAL;
      else if( "STEREOPARENT".equalsIgnoreCase(molOpt) )
         outMol = OUTMolFormat.STEREOPARENT;
      else if( "NORMALIZED".equalsIgnoreCase(molOpt) )
         outMol = OUTMolFormat.NORMALIZED;
      else if( "TAUTOMERIC".equalsIgnoreCase(molOpt) )
         outMol = OUTMolFormat.TAUTOMERIC;
      else {
         System.err.printf("Unkown option for -mol: %s\n", molOpt);
         System.exit(1);
      }

      String inFile  = cmd.getOptionValue("in");
      String outFile = cmd.getOptionValue("out");
      boolean limitMessage = cmd.hasOption("shortMessage");

      try {
         oemolistream ifs = new oemolistream(inFile);
         oemolostream ofs = new oemolostream(outFile);

         URL cFile;
         if( cmd.hasOption("f") )
         {  File cf = new File(cmd.getOptionValue("f"));
            if( ! cf.canRead() ) throw new Error("Can not read from: " + cf.getAbsolutePath());
               cFile= cf.toURI().toURL();
         } else
         {  cFile = OEStruchk.getResourceURL(OEStruchk.class,"Struchk.xml");
         }

         // create OEStruchk from config file
         OEStruchk strchk = new OEStruchk(cFile, CHECKConfig.ASSIGNStructFlag, false);

         OEGraphMol mol   = new OEGraphMol();
         StringBuilder sb = new StringBuilder(2000);
         while ( oechem.OEReadMolecule(ifs , mol ) ) {
            int nNonTetrahedralChiral = 0;
            String d = oechem.OEGetSDData(mol, NON_TETRAHEDRAL_CHIRAL_TAG);
            if( d != null && d.length() > 0 ) nNonTetrahedralChiral = Integer.parseInt(d);

            if(! strchk.applyRules(mol, null,nNonTetrahedralChiral))
               nErrors++;

            switch( outMol ) {
               case STEREOPARENT:
                  mol.Clear();
                  oechem.OEAddMols(mol, strchk.getTransformedMol("parentAllStereo"));
                  break;
               case NORMALIZED:
                  mol.Clear();
                  oechem.OEAddMols(mol, strchk.getTransformedMol("parent"));
                  break;
               case TAUTOMERIC:
                  mol.Clear();
                  oechem.OEAddMols(mol, strchk.getTransformedMol(null));
                  break;
               case ORIGINAL:
               break;
            }

            oechem.OESetSDData(mol, "CTISMILES", strchk.getTransformedIsoSmiles(null));
            oechem.OESetSDData(mol, "CTSMILES",  strchk.getTransformedSmiles(null));
            oechem.OESetSDData(mol, "CISMILES",  strchk.getTransformedIsoSmiles("parent"));
            oechem.OESetSDData(mol, "CTSSMILES", strchk.getTransformedIsoSmiles(
                                                   StructFlagAnalysisInterface.STEREONormalizedKeeper));
            oechem.OESetSDData(mol, "Strutct_Flag",  strchk.getStructureFlag().getName());

            List<Message> msgs = strchk.getStructureMessages(null);
            nMessages += msgs.size();
            for( Message msg : msgs )
               sb.append(String.format("\t%s:%s", msg.getLevel(), msg.getText()));
            if( limitMessage ) sb.setLength( Math.min(sb.length(), 80));

            oechem.OESetSDData(mol, "NORM_MESSAGE",  sb.toString());
            oechem.OEMDLPerceiveBondStereo(mol);
            oechem.OEWriteMolecule(ofs, mol);

            sb.setLength(0);
            nStruct++;
         }
         strchk.delete();
         mol.delete();
         ifs.close();
         ifs.delete();
         ofs.close();
         ofs.delete();

      } catch (Exception e) {
         throw new Error(e);
      }finally {
         System.err.printf("sdfNormalizer: Checked %d structures %d errors, %d messages in %dsec\n",
               nStruct, nErrors, nMessages, (System.currentTimeMillis()-start)/1000);
      }
   }

   private static void exitWithHelp(Options options, String msg)
   {  System.err.println(msg);
      HelpFormatter formatter = new HelpFormatter();
      formatter.printHelp("sdfNormalizer", options);
      System.err.println();
      System.err.println("Will validate that the structures int he input file conform to the \n"
           + "Genentech normalization rules. The normalization also includes stripping of known salts and solvents.\n"
           + "The output file will include the isomeric and non-simeric smiles\n"
           + "For the tautomer of the input and for the unique tautomer generated by quackpac\n"
           + "The -mol option the molfile of each record will contain the unmodified input or the normalized structture\n");
      System.err.println();

      System.exit(1);
   }


}
