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

import openeye.oechem.OEAtomBase;
import openeye.oechem.OEGraphMol;
import openeye.oechem.OEIsHeavy;
import openeye.oechem.oechem;
import openeye.oechem.oemolistream;
import openeye.oechem.oemolostream;

import org.apache.commons.cli.*;

public class SDFFilter
{
   private static int maxAtoms = 100;

//   private static boolean passFilters(OEGraphMol mol, OEFilter filter)
   private static boolean passFilters(OEGraphMol mol)
   {
      //number of heavy atoms should be less than maxAtoms                                                                                                                        
      if (oechem.OECount(mol, new OEIsHeavy()) > maxAtoms) {
         System.err.println("Ignoring molecule name: " + mol.GetTitle() + ", it has more than " + maxAtoms + " atoms.");
         return false;
      }
      //number of heavy atoms should be greater than 1
      if (oechem.OECount(mol, new OEIsHeavy()) < 1) {
         System.err.println("Ignoring molecule name: " + mol.GetTitle() + ", it does not have any atoms");
         return false;
      }
      //keep molecules containing one component
      int[] parts = new int[mol.GetMaxAtomIdx()];
      int count = oechem.OEDetermineComponents(mol, parts);
      if (count > 1) {
         System.err.println("Ignoring molecule name: " + mol.GetTitle() + ", it has " + count + " structures");
         return false;
      }

      //atomic elements greater than 1
      for (OEAtomBase atom : mol.GetAtoms()) {
         if (atom.GetAtomicNum() < 1) {
            System.err.println("Ignoring molecule name: " + mol.GetTitle() + ", it has an invalid atom: " + atom.GetName());
            return false;
         }
         //ignore atoms greater than Iodine
         if (atom.GetAtomicNum() > 53) {
            System.err.println("Ignoring molecule name: " + mol.GetTitle() + ", it has an atomic number greater than 53 (Iodine): " + atom.GetName());
            return false;
         }
      }

      //filter out compounds with invalid valence
      /* it appears that we don't have license to the filter toolkit
      if (! filter.call(mol)) {
         System.err.println("Failed valence check: " + mol.GetTitle() );
         return false;
      }
      */

      return true;
   }

   public static void main(String args[]){
      String usage = "java sdfFilter [options]\n" +
      		         "Filters out molecules that contain the following:\n" + 
                     "1. More than maxHeavy atoms (default is 100)\n" +
      		         "2. Have more than one component\n"+
                     "3. Does not have any atoms\n"+
      		         "4. Have invalid atoms, like R, * etc ...\n" +
                     "5. Have atomic number greater than 53 (Iodine)";

      // create Options object
      Options options = new Options();
      // add  options
      options.addOption("h", false, "");
      options.addOption("in", true, "inFile in OE formats: Ex: a.sdf or .sdf");
      options.addOption("out", true, "outFile in OE formats. Ex: a.sdf or .sdf");
      options.addOption("filtered", true, "Optional: filteredFile in OE formats. Ex: a.sdf or .sdf");
      options.addOption("maxHeavy", true, "Number of heavy atom cutoff. Default is 100");
      CommandLineParser parser = new PosixParser();

      try
      {
         CommandLine cmd = parser.parse(options, args);
         String inFile = cmd.getOptionValue("in");
         if (inFile == null) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp( usage, options );
            System.exit(1);
         }

         String outFile = cmd.getOptionValue("out");
         if (outFile == null) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp( usage, options );
            System.exit(1);
         }

         String filteredFile=null;
         if (cmd.hasOption("filtered")) {
            filteredFile = cmd.getOptionValue("filtered");
         }
         
         int maxHeavy = 100;
         if (cmd.hasOption("maxHeavy")) {
            maxHeavy = Integer.parseInt(cmd.getOptionValue("maxHeavy"));
            maxAtoms=maxHeavy;
         }

         if (cmd.hasOption("h")) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp( usage, options );
            System.exit(1);
         }
         
         oemolistream ifs = new oemolistream(inFile);  
         oemolostream ofs = new oemolostream(outFile);  
         oemolostream ffs = null;
         if (filteredFile != null) {
            ffs= new oemolostream(filteredFile);  
         }
         /*it appears that we don't have license to molprop toolkit
         oeisstream iss = new oeisstream();
         OEFilter filter = new OEFilter(iss);
         filter.SetTypeCheck(true);
         oechem.OEThrow.SetLevel(OEErrorLevel.Warning);
         */

         OEGraphMol mol = new OEGraphMol();

         while (oechem.OEReadMolecule(ifs, mol))
         {
//            if (passFilters(mol, filter)) { //needs license to molprop toolkit
            if (passFilters(mol)) {
               oechem.OEWriteMolecule(ofs, mol);
            }else {
               if (ffs != null) {
                  oechem.OEWriteMolecule(ffs, mol);
               }
            }
         }

      } catch (ParseException e)
      {  // TODO print explaination
         throw new Error(e);
      } 
   }
}
