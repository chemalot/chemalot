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
package com.genentech.chemistry.tool;

import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

import openeye.oechem.*;
import openeye.oedepict.OEAlignmentOptions;
import openeye.oedepict.OEAlignmentResult;
import openeye.oedepict.oedepict;

import org.apache.commons.cli.*;

import com.genentech.oechem.tools.OETools;

/**
 * @author jianwenf
 *
 */
public class sdf2DAlign
{
   private static final String MY_NAME = "sdf2DAlign \n" +
               "Aligns 2D coordinates of input file to substructures specified in " +
               "the templates file. Once a substructure is found, the rest of the " +
               "substructure in the templates file are ignored.";

   public static final OEAlignmentOptions ALIGNOpts = new OEAlignmentOptions();

   //class variables
   private String inFile, outFile, templateFile;
   private boolean debug;

   public void setDebug(boolean debug)
   {
      this.debug = debug;
   }

   //constructor
   /**
    * align depictions of molecules in inFile to templates in templateFile
    * Results are written out to outFile
    * resetCoords provide option to reset 2D coordinates using OE algorithm
    * @param inFile
    * @param outFile
    * @param templateFile
    * @param resetCoords
    */
   public sdf2DAlign (String inFile, String outFile, String templateFile, boolean resetCoords)
   {
      this.inFile=inFile;
      this.outFile=outFile;
      this.templateFile=templateFile;
      setAlignmentOptions(resetCoords);
   }

   private static void exitWithHelp(Options options)
   {
      HelpFormatter formatter = new HelpFormatter();
      formatter.printHelp(MY_NAME, options);
      System.exit(1);
   }

   /**
    * Read templates file and returns a list of Substructure search objects
    * Template file is set by constructor method
    * @return
    */
   private List<OESubSearch> getTemplates()
   {
      oemolistream ifs = new oemolistream();
      ifs.SetFormat(OEFormat.MDL);
      int aromodel = OEIFlavor.Generic.OEAroModelOpenEye;
      int qflavor  = ifs.GetFlavor( ifs.GetFormat());
      ifs.SetFlavor(ifs.GetFormat(),(qflavor|aromodel));

      if (!ifs.open(this.templateFile))
         oechem.OEThrow.Fatal("Unable to create " + templateFile);

      List <OESubSearch> templates = new LinkedList<OESubSearch>();
      int opts = OEMDLQueryOpts.Default|OEMDLQueryOpts.SuppressExplicitH;
      OEQMol qmol = new OEQMol();
      while (oechem.OEReadMDLQueryFile(ifs, qmol, opts))
      {
         System.err.println("Reading template molecule " + qmol.GetTitle() + " with SMILES: " + oechem.OECreateSmiString(qmol));
         qmol.SetDimension(2);
         OESubSearch ss = new OESubSearch(qmol);
//         ss.SetMaxMatches(1);

         templates.add(ss);
         qmol.Clear();
      }
      qmol.delete();
      ifs.close();
      ifs.delete();
      return templates;
   }

   /**
    * changes molecules in place
    * @param templates
    */
   private void applyTemplates(List<OESubSearch> templates)
   {
      oemolistream ifs = new oemolistream();
      oemolostream ofs = new oemolostream();

      if (!ifs.open(this.inFile))
         oechem.OEThrow.Fatal("Unable to open " + this.inFile);
      if (!ofs.open(this.outFile))
         oechem.OEThrow.Fatal("Unable to create " + this.outFile);

      OEGraphMol mol = new OEGraphMol();
      while (oechem.OEReadMolecule(ifs, mol))
      {
        alignDepictionToFirstTemplate(templates, mol);
        oechem.OEWriteMolecule(ofs, mol);
        mol.Clear();
      }
      mol.delete();
      ifs.close();
      ofs.close();
      ifs.delete();
      ofs.delete();
   }

   /**
    * Align depictions to first template that matches
    * @param templates
    * @param mol
    */
   private void alignDepictionToFirstTemplate(List<OESubSearch> templates, OEGraphMol mol)
   {  // The next line (OEPrepareDepiction) should not be required, but as of 2016.Feb without this oedepict will crash
      oedepict.OEPrepareDepiction(mol, ALIGNOpts.GetClearCoords(), true);

      //loop through list of substructures and stop on first match
      //option to reset coordinates if specified in command line
      for (OESubSearch ss : templates)
      {  OEAlignmentResult aRes = oedepict.OEPrepareAlignedDepiction(mol, ss, ALIGNOpts);
         boolean matches = aRes.IsValid();
         if (matches && this.debug)
            System.err.printf("%s matched %s with %frmsd\n",
                     OETools.molToCanSmi(mol,true),ss.GetPattern().GetTitle(), aRes.GetRMSD());
         aRes.delete();
         if( matches )
         {  oechem.OEMDLPerceiveBondStereo(mol); // requirted to fix stereochemsitry bug that should be solved in June 2017
            return;
         }
      }
      //did not find any substructure matches if code gets here
      if( debug )
         System.err.println("Failed to find a substructure match for " + mol.GetTitle() + " :" + oechem.OECreateSmiString(mol));
   }


   private static void setAlignmentOptions(boolean resetCoords)
   {  ALIGNOpts.SetClearCoords(resetCoords);
      ALIGNOpts.SetRotateAroundBonds(true);
      ALIGNOpts.SetMaxBondRotations(3200);
      ALIGNOpts.SetFixedCoords(true);
      ALIGNOpts.SetRMSDCutoff(0.1);
      ALIGNOpts.SetMaxPatternMatches(20);
   }



   public static void main(String[] args) throws IOException
   {

      final String OPT_INFILE         = "in";
      final String OPT_OUTFILE        = "out";
      final String OPT_TEMPLATEFILE   = "templates";
      final String OPT_RESETCOORDS    = "reset";
      final String OPT_DEBUG          = "debug";

      Options options = new Options();

      Option opt = new Option(OPT_INFILE, true,
               "Input file (oe-supported). Use .sdf to specify the file type.");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option(OPT_OUTFILE, true,
               "Output file (oe-supported). Use .sdf to specify the file type.");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option(OPT_TEMPLATEFILE, true,
               "Template file (oe-supported). Use .sdf to specify the file type.");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option(OPT_RESETCOORDS, false,
               "Reset coordinates in input file when aligning. This will delete original 2D coords and create new ones.");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option(OPT_DEBUG, false,
               "Print debug statements.");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("h", false,
               "print usage statement");
      opt.setRequired(false);
      options.addOption(opt);

      PosixParser parser = new PosixParser();
      CommandLine cl=null;

      try {
         cl = parser.parse(options, args);
      } catch (Exception exp) {
         System.err.println(exp.getMessage());
         exitWithHelp(options);
      }

      //get command line parameters
      String inFile = cl.getOptionValue(OPT_INFILE);
      String outFile = cl.getOptionValue(OPT_OUTFILE);
      String templateFile =cl.getOptionValue(OPT_TEMPLATEFILE);
      boolean resetCoords =cl.hasOption(OPT_RESETCOORDS);
      boolean debug =cl.hasOption(OPT_DEBUG);

      sdf2DAlign myAlign = new sdf2DAlign(inFile, outFile, templateFile, resetCoords );
      myAlign.setDebug(debug);

      //get templates as subsearches
      List<OESubSearch> templates = myAlign.getTemplates();

      //apply template to input molecules
      myAlign.applyTemplates(templates);

      for(OESubSearch ss : templates) ss.delete();
   }
}

