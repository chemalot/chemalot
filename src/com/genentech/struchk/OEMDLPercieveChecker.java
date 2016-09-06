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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

import openeye.oechem.*;

import org.apache.commons.cli.*;
import org.jdom.JDOMException;

import com.genentech.oechem.tools.OETools;

/**
 * SlatsTripper removes the salts and adds the saltCode MW and MF to an sdf file.
 * @author albertgo 2008
 *
 */
public class OEMDLPercieveChecker {

   private OEGraphMol tmpMol;

   OEMDLPercieveChecker( ) {
      tmpMol = new OEGraphMol();
   }
   
   /**
    * Strip the mol and return a stripped mol object with associated data in the
    * SD tags.
    * 
    * @param mol input mol.
    * @return 
    * @return transformed molecule object, do not modify or delete it is owned by SaltStripper.
    */
   boolean checkMol(OEGraphMol mol) {
      tmpMol.Clear();
      oechem.OEMDLStereoFromBondStereo(mol);
      oechem.OEAddMols(tmpMol, mol);
      oechem.OEMDLPerceiveBondStereo(mol);
      oechem.OEMDLCorrectBondStereo(tmpMol);
      
      String molSmi = OETools.molToCanSmi(mol, true);
      String corrSmi = OETools.molToCanSmi(tmpMol, true);
      if(molSmi.equals(corrSmi)) return true;
      
      System.err.printf("%s>>%s\n", molSmi, corrSmi);
      oechem.OESetSDData(mol, "FromBondsSmi", molSmi);
      oechem.OESetSDData(mol, "FromParitySmi", corrSmi);
      
      return false;
   }
   
   /** after calling delete this SaltStripper should not be used again */
   public void delete() {
      tmpMol.delete();
   }
   
   public static void main(String [] args) throws ParseException, JDOMException, IOException {
      // create command line Options object
      Options options = new Options();
      Option opt = new Option("i",true, "input file [.ism,.sdf,...]");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("o",true, "output file");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("d",false, "debug: wait for user to press key at startup.");
      opt.setRequired(false);
      options.addOption(opt);

      CommandLineParser parser = new PosixParser();
      CommandLine cmd = null;
      try {
         cmd = parser.parse( options, args);
      } catch(Exception e)
      {  System.err.println(e.getMessage());
         exitWithHelp(options);
      }
      args = cmd.getArgs();

      if(cmd.hasOption("d")) {
         System.err.println("Start debugger and press return:");
         new BufferedReader(new InputStreamReader(System.in)).readLine();
      }
       
      if(args.length != 0) { 
         exitWithHelp(options);
      }
      
      String inFile  = cmd.getOptionValue("i");
      String outFile = cmd.getOptionValue("o");
      
      OEMDLPercieveChecker checker = null;
      try {
         checker = new OEMDLPercieveChecker( );
         
         oemolostream out = new oemolostream(outFile);
         oemolistream in = new oemolistream(inFile);
         
         OEGraphMol mol = new OEGraphMol();
         while ( oechem.OEReadMolecule(in , mol ) ) {
            if( ! checker.checkMol(mol) )
            oechem.OEWriteMolecule(out, mol);
         }
         checker.delete();
         in.close();
         in.delete();
         
         out.close();
         out.delete();
         
      } catch (Exception e) {
         throw new Error(e);
      }
      System.err.println("Done:");
   }

   private static void exitWithHelp(Options options) {
      HelpFormatter formatter = new HelpFormatter();
      formatter.printHelp( "OEMDLPercieveChecker", options );
      System.exit(1);
   }
}
