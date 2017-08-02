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

package com.genentech.chemistry.tool.align;

import java.util.ArrayList;

import openeye.oechem.*;

import org.apache.commons.cli.*;

public class SDFAlign
{  public enum OUTType
   {  ALIGNED,ORIGINAL; }

   public static boolean is3DFormat(int fmt)
   {  if (fmt == OEFormat.SMI || fmt == OEFormat.ISM || fmt == OEFormat.CAN|| fmt == OEFormat.MF)
      {  return false;
      }
      return true;
   }

   public static void main(String...args)
   {  Options options = new Options();
      Option opt = new Option("in", true, "input sd file");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("out", true, "output file");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("method", true, "fss|sss|MCS|clique (default mcs).");
      options.addOption(opt);

      opt = new Option("ref", true, "reference molecule (if not given first molecule in '-in' is used). If multiple ref molecules are read the min RMSD is reported.");
      options.addOption(opt);

      opt = new Option("mirror", false, "If given and the molecule is not chiral, return best mirror image.");
      options.addOption(opt);

      opt = new Option("rmsdTag", true, "Tagname for output of rmsd, default: no output.");
      options.addOption(opt);

      opt = new Option("atomMatch", true, "Sequence of none|default|hcount|noAromatic specifing how atoms are matched cf. oe document.\n"
                                         +"noAromatic can be used to make terminal atoms match aliphatic and aromatic atoms.\n"
                                         +"Queryfeatures are considered only if default is used.");
      options.addOption(opt);

      opt = new Option("bondMatch", true, "Sequence of none|default specifing how bonds are matched cf. oe document.");
      options.addOption(opt);

      opt = new Option("keepCoreHydrogens", false, "If not specified the hydrigen atoms are removed from the core.");
      options.addOption(opt);

      opt = new Option("outputMol", true, "aligned|original (def: aligned) use original to just compute rmsd.");
      options.addOption(opt);

      opt = new Option("doNotOptimize", false, "If specified the RMSD is computed without moving optimizing the overlay.");
      options.addOption(opt);

      opt = new Option("quiet", false, "Reduced warining messages");
      options.addOption(opt);

      CommandLineParser parser = new BasicParser();
      CommandLine cmd = null;
      try
      {
         cmd = parser.parse(options, args);
      } catch (Exception e)
      {  exitWithHelp(e.getMessage(), options);
      }

      if (cmd.getArgs().length > 0)
         exitWithHelp("To many arguments", options);

      // do not check aromaticity on atoms so that a terminal atom matches aromatic and non aromatic atoms
      int atomExpr = OEExprOpts.DefaultAtoms;
      int bondExpr = OEExprOpts.DefaultBonds;

      String atomMatch = cmd.getOptionValue("atomMatch");
      if( atomMatch == null ) atomMatch = "";
      atomMatch = '|' + atomMatch.toLowerCase() + '|';

      String bondMatch = cmd.getOptionValue("bondMatch");
      if( bondMatch == null ) bondMatch = "";
      bondMatch = '|' + bondMatch.toLowerCase() + '|';

      String inFile = cmd.getOptionValue("in");
      String outFile = cmd.getOptionValue("out");
      String refFile = cmd.getOptionValue("ref");
      String method = cmd.getOptionValue("method");
      String rmsdTag = cmd.getOptionValue("rmsdTag");
      String oMol = cmd.getOptionValue("outputMol");
      boolean doMirror = cmd.hasOption("mirror");
      boolean doOptimize = !cmd.hasOption("doNotOptimize");
      boolean quiet = cmd.hasOption("quiet");

      OUTType outputMol = oMol == null ? OUTType.ALIGNED : OUTType.valueOf(oMol.toUpperCase());

      if( atomMatch.startsWith("|none") )      atomExpr = 0;
      if( atomMatch.contains("|hcount|") )     atomExpr |= OEExprOpts.HCount;
      if( atomMatch.contains("|noAromatic|") ) atomExpr &= (~ OEExprOpts.Aromaticity);

      if( bondMatch.startsWith("|none") ) bondExpr = 0;

      ArrayList<OEMol> refmols = new ArrayList<OEMol>();
      if( refFile != null )
      {  oemolistream reffs = new oemolistream(refFile);
         if (!is3DFormat(reffs.GetFormat()))
            oechem.OEThrow.Fatal("Invalid input format: need 3D coordinates");
         reffs.SetFormat(OEFormat.MDL);

         int aromodel = OEIFlavor.Generic.OEAroModelOpenEye;
         int qflavor  = reffs.GetFlavor( reffs.GetFormat());
         reffs.SetFlavor(reffs.GetFormat(),(qflavor|aromodel));

         OEMol rmol = new OEMol();
         while(oechem.OEReadMDLQueryFile(reffs,rmol))
         {  if( ! cmd.hasOption("keepCoreHydrogens"))
               oechem.OESuppressHydrogens(rmol);
            refmols.add(rmol);
            rmol = new OEMol();
         }
         rmol.delete();

         if( refmols.size() == 0)
            throw new Error("reference file had no entries");

         reffs.close();
      }

      oemolistream fitfs = new oemolistream(inFile);
      if (!is3DFormat(fitfs.GetFormat()))
         oechem.OEThrow.Fatal("Invalid input format: need 3D coordinates");

      oemolostream ofs = new oemolostream(outFile);
      if (!is3DFormat(ofs.GetFormat()))
         oechem.OEThrow.Fatal("Invalid output format: need 3D coordinates");

      AlignInterface aligner = null;
      OEGraphMol fitmol = new OEGraphMol();

      if( oechem.OEReadMolecule(fitfs, fitmol) )
      {  if( refmols.size() == 0)
         {  OEMol rmol = new OEMol(fitmol);
            if( ! cmd.hasOption("keepCoreHydrogens"))
               oechem.OESuppressHydrogens(rmol);

            refmols.add(rmol);
         }

         if( "sss".equalsIgnoreCase(method))
         {  aligner = new SSSAlign(refmols, outputMol, rmsdTag, doOptimize, doMirror, atomExpr, bondExpr, quiet);

         } else if( "clique".equalsIgnoreCase(method))
         {  aligner = new CliqueAlign(refmols, outputMol, rmsdTag, doOptimize, doMirror, atomExpr, bondExpr, quiet);

         }else if( "fss".equalsIgnoreCase(method))
         {  if( cmd.hasOption("atomMatch") || cmd.hasOption("bondMatch") )
               exitWithHelp("method fss does not support '-atomMatch' or '-bondMatch'", options);
            aligner = new FSSAlign(refmols, outputMol, rmsdTag, doOptimize, doMirror);

         } else
         {  aligner = new McsAlign(refmols, outputMol, rmsdTag, doOptimize, doMirror, atomExpr, bondExpr, quiet);
         }

         do
         {  aligner.align(fitmol);
            oechem.OEWriteMolecule(ofs, fitmol);
         }while (oechem.OEReadMolecule(fitfs, fitmol));

      }

      fitmol.delete();
      if( aligner != null ) aligner.close();
      for( OEMolBase mol : refmols) mol.delete();
      fitfs.close();
      ofs.close();
   }

   private static void exitWithHelp(String msg, Options options)
   {  System.err.println(msg);

      HelpFormatter formatter = new HelpFormatter();
      formatter.setWidth(120);
      formatter.printHelp( "sdfAlign -in .sdf -out .sdf [-ref .sdf]", options );
      System.exit(1);
   }

}
