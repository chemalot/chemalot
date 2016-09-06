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

/**
 * This class provides the framework to call aggregation functions such as count, mean, etc.
 * @author Johnny Wu Aug 11, 2011
 *
 */

import java.io.IOException;
import java.util.regex.Pattern;

import openeye.oechem.OEGraphMol;
import openeye.oechem.oechem;
import openeye.oechem.oemolistream;
import openeye.oechem.oemolostream;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.PosixParser;

public class SDFBinning
{
   private static String inFile;
   private static String outFile;
   private static double[] binValues;
   private static String[] binNames;
   private static String binTag;
   private static String dataTag;
   private static oemolostream ofs;
   private static oemolistream ifs;
   private static boolean ignoreModifier;

   private static final Pattern REMOVEModPat = Pattern.compile("^[<> ~=]+");


   private static void exitWithHelp(Options options)
   {
      HelpFormatter formatter = new HelpFormatter();
      formatter.printHelp("SDFBinning", options);
      System.exit(1);
   }


   private static void processOptions(String[] args)
   {
      Options options = new Options();

      Option opt = new Option("in", true,
               "Input file oe-supported Use .sdf to specify the file type.");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("out", true,
               "Output file oe-supported Use .sdf to specify the file type.");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("bin", true,
               "bin definition. Use format \"bin_1_upperbound=bin_1_name|bin_2_upperbound=bin_2_name|bin_3_upperbound_=bin_3_name\"");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("dataTag", true,
               "SD tag containing data for binning");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("binTag", true,
               "SD tag used for storing bin names");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("ignoreModifier", false,
               "ignore leading '><= ~'");
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
         // catch (ParseException exp) {
         System.err.println(exp.getMessage());
         exitWithHelp(options);
      }
      inFile = cl.getOptionValue("in");
      outFile = cl.getOptionValue("out");
      ignoreModifier = cl.hasOption("ignoreModifier");
      String binDefinition=cl.getOptionValue("bin");
      parseBins(binDefinition);

      binTag=cl.getOptionValue("binTag");
      dataTag=cl.getOptionValue("dataTag");
   }

   //-bin "6.2=S|15=M|100000=L" -binTag "HLM CLhep [ml/min/kg]"
   private static void parseBins(String binDefinition)
   {// split by vertical |, set binValues and binNames
      String[] bins = binDefinition.split("\\|");
      binValues = new double[bins.length];
      binNames = new String[bins.length];


      try {
         for (int i=0; i<bins.length; i++){
            String [] valueNamePair = bins[i].split("=");
            if (valueNamePair.length == 2){
               binValues[i]=Double.parseDouble(valueNamePair[0]);
               binNames[i]=valueNamePair[1];
            }else {
               throw new Exception();
            }
         }
      } catch (Exception e) {
         System.err.println("Invalid bin definition: " + binDefinition);
         System.exit(1);
      }
   }

   /*
    * Process each molecule
    */
   private static void processRecord(OEGraphMol mol)
   {
      // Process current record
      // Get data from dataTag, assign bin, save bin name in binTag
      String value = oechem.OEGetSDData(mol, dataTag);
      String binName = assignBin(value);
      oechem.OESetSDData(mol, binTag, binName);
      oechem.OEWriteMolecule(ofs, mol);
   }
/*returns "" if a bin cannot be assigned */
   private static String assignBin(String value)
   {
      try{
         if( ignoreModifier ) value = REMOVEModPat.matcher(value).replaceFirst("");

         double val = Double.parseDouble(value);
         for (int i=0; i<binValues.length; i++){
            if (val < binValues[i]){
               return binNames[i];
            }
         }
      }catch (Exception e){
         return "";
      }
      return "";
   }

   public static void main(String[] args) throws IOException
   {

      processOptions(args);

      ifs = new oemolistream();
      ofs = new oemolostream();



      if (!ifs.open(inFile))
         oechem.OEThrow.Fatal("Unable to open " + inFile);
      if (!ofs.open(outFile))
         oechem.OEThrow.Fatal("Unable to create " + outFile);

      OEGraphMol mol = new OEGraphMol();
      while (oechem.OEReadMolecule(ifs, mol))
      {
         processRecord(mol);
      }

      ifs.close();
      ofs.close();
   }
}
