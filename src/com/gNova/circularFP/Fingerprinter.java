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

package com.gNova.circularFP;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import static com.gNova.circularFP.CFingerprint.CFPCountType;

import openeye.oechem.*;

import org.apache.commons.cli.*;
/**
 *
 * @author tjodonnell
 *
 */

public class Fingerprinter
{
   private static String infile = null;
   private static String outfile = null;
   private static String format = null;
   private static String type = null;
   private static String smaFile = null;

   //private statis String aromaticity = null;
   private static int[] levels = new int[] {0,1,2,3};
   private static int nbits = 256;
   private static boolean verbose = false;
   private static CFPCountType countType;

   public static void usage(String msg, Options options)
   {  if( msg != null && msg.length()>0)
         System.err.println(msg);
      HelpFormatter helper = new HelpFormatter();
      helper.printHelp( "sdfCFP", options );
      System.exit(0);
   }

   private static Options defineCommandLineOptions()
   {
      // create Options object
      Options options = new Options();

      // add t option
      options.addOption("help",    false, "print this message");
      options.addOption("in",      true,  "REQUIRED input file.  Use .sdf (.smi, etc.) for stdin");
      options.addOption("out",     true,  "output file: default .sdf (to stdout)");
      options.addOption("level",   true,  "fingerprint radius, repeats allowed. default 0-3");
      options.addOption("format",  true,  "output format: hex, bitlist, atomID, counts: default hex");
      options.addOption("functionDefinition",  true, "smarts file defining inital functional atom identifiers");
      //options.addOption("aromaticity",  true, "aromaticity model: oechem, mdl: default oechem");
      options.addOption("nbits",   true,  "fingerprint length: default 256");
      options.addOption("count",   true,  "For atomId, bitlist and hex format: set additional bits for counts options: NOCount|LOGCount|LINCount");
      options.addOption("type",    true,  "fingerprint type: atomic, functional: default atomic");
      options.addOption("verbose", false, "verbose output, for debugging");
      return options;
   }

   private static void validateOptions(Options options)
   {
      if (!(format.equals("hex")
            || format.equals("counts")
            || format.equals("bitlist")
            || format.equals("atomID")))
         usage("Invalid format: " + format, options);

      if( countType != CFPCountType.NOCount && !(format.equals("hex") || format.equals("bitlist")))
         usage("count only accpeted or hex or bitlist format", options);

      if (!(type.equals("atomic")
          || type.equals("functional")))
         usage("Invalid type: " +type, options);
   }

   private static void parseCommandLine(String...args)
   {
      CommandLineParser parser = new PosixParser();
      Options options = defineCommandLineOptions();
      try
      {  CommandLine cmd = parser.parse( options, args);
         if (cmd.hasOption("help")) usage(null, options);
         if (!cmd.hasOption("in")) usage("In file required", options);

         infile = cmd.getOptionValue("in", ".sdf");
         outfile = cmd.getOptionValue("out", ".sdf");
         format = cmd.getOptionValue("format", "hex");
         type = cmd.getOptionValue("type", "atomic");
         smaFile = cmd.getOptionValue("functionDefinition");

         countType = CFPCountType.NOCount;
         if(cmd.hasOption("count"))
            countType = CFPCountType.valueOf(cmd.getOptionValue("count"));

         if (cmd.hasOption("verbose")) Fingerprinter.verbose = true;

         String[] lvlStr = cmd.getOptionValues("level");
         if( lvlStr != null )
         {  levels = new int[lvlStr.length];
            for(int i=0; i< lvlStr.length; i++)
               levels[i] = Integer.parseInt(lvlStr[i]);
         }

         nbits = Integer.parseInt(cmd.getOptionValue("nbits", "256"));
      } catch (Exception exp)
      {
         String msg = "Parsing failed: " + exp.getMessage();
         usage(msg, options);
      }
      validateOptions(options);
   }


   public static void main(String...args)
   {  parseCommandLine(args);

      // option to fp.generate for verbose output

      CFP ecfp = new CFP(verbose);
      oemolistream ifs = new oemolistream(Fingerprinter.infile);
      oemolostream ofs = new oemolostream(Fingerprinter.outfile);
      OEGraphMol mol = new OEGraphMol();

   /*
      // read first molecule; Tanimoto of others uses this initial fp
      oechem.OEReadMolecule(ifs, mol);
      ecfp.generate(mol, 3, verbose);
      Fingerprint fp0 = new CFingerprint(ecfp.get(), 512);
      //System.out.print(String.format("%5.3f",fp0.Tanimoto(fp0)));
      System.out.print(fp0.getNBits()+" ");
      oechem.OEWriteMolecule(ofs, mol);
      mol.Clear();
   */

      String tag;
       if (type.equals("functional"))
      {  tag = "FFP";
         List<String> smarts = null;

         if (smaFile != null)
            smarts = readSmartsFile();
         ecfp.initializeSmarts(smarts);

      } else
      {  tag = "AFP";
      }

      while (oechem.OEReadMolecule(ifs, mol))
      {  Fingerprint fp;

         for (int lvl : levels )
         {  ecfp.clear();
            ecfp.generate(mol, lvl, type);
            if (format.equals("hex"))
            {  fp = CFingerprint.createCFingerprint(ecfp.getCounts(0), nbits, countType);
               oechem.OESetSDData(mol, tag+lvl, fp.getHexString());

            } else if (format.equals("bitlist"))
            {  fp = CFingerprint.createCFingerprint(ecfp.getCounts(0), nbits, countType);
               oechem.OESetSDData(mol, tag+lvl, fp.getBitString());

            } else if (format.equals("counts"))
            {  Map<Integer,Integer> counts = ecfp.getCounts(nbits);
               String myTag = tag+lvl + '_';
               for(int i=0; i< nbits; i++)
               {  Integer count = counts.get(i);
                  if( count != null )
                     oechem.OESetSDData(mol, (myTag + i), counts.get(i).toString());
                  else
                     oechem.OESetSDData(mol, (myTag + i), "0");
               }
            } else if (format.equals("atomID"))
            {  fp = CFingerprint.createCFingerprint(ecfp.getCounts(0), nbits, countType);
               oechem.OESetSDData(mol, tag+lvl, fp.getAtomIDString());

            }
            //tag = "CFPBits" + i;
            //oechem.OESetSDData(mol, tag, Integer.toString(fp.getNBits()));
         }
         oechem.OEWriteMolecule(ofs, mol);
         mol.Clear();
      }
      ifs.close();
      ofs.close();
   }

   private static List<String> readSmartsFile()
   {  ArrayList<String> smarts = new ArrayList<String>();
      try
      {  FileInputStream fstream = new FileInputStream(smaFile);
         DataInputStream in = new DataInputStream(fstream);
         BufferedReader br = new BufferedReader(new InputStreamReader(in));
         String strLine;
         while ((strLine = br.readLine()) != null)
         {  if (!strLine.startsWith("#") && strLine.trim().length() > 0)
            {
               smarts.add(strLine);
            }
         }
         in.close();
      } catch (Exception e){
         System.err.println("Error: " + e.getMessage());
         System.exit(0);
      }
      return smarts;
   }
}
