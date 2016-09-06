/*
   Copyright 2006-2014 Man-Ling Lee & Alberto Gobbi

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

Contact: aestelSW@gmail.com
*/

package com.aestel.chemistry.openEye.fp;

import java.io.*;
import java.util.Random;

import openeye.oechem.*;

import org.apache.commons.cli.*;

/**
 * Sorts the indexes of the fragments in the dictionary so that frequent fragments
 * have smaller indexes. This is to save space and speed up oracle search.
 *  @author albertgo
 *
 */
public class FPDictionarySorter
{  public static void main(String...args) throws IOException
   {  long start = System.currentTimeMillis();
      int iCounter = 0;
      int fpCounter = 0;
   
      // create command line Options object
      Options options = new Options();
      Option opt = new Option("i",true, "input file [.ism,.sdf,...]");
      opt.setRequired(true);
      options.addOption(opt);
   
      opt = new Option("fpType",true, "fingerPrintType: maccs|linear7|linear7*4");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("sampleFract",true, 
                        "fraction of input molecules to use (Default=1)");
      opt.setRequired(false);
      options.addOption(opt);

      CommandLineParser parser = new PosixParser();
      CommandLine cmd = null;
      try
      {  cmd = parser.parse( options, args);
      } catch(Exception e)
      {  System.err.println(e.getMessage());
         exitWithHelp(options);
      }
      args = cmd.getArgs();

      if(cmd.hasOption("d"))
      {  System.err.println("Start debugger and press return:");
         new BufferedReader(new InputStreamReader(System.in)).readLine();
      }

      if(args.length != 0)
      {  exitWithHelp(options);
      }
      
      String type = cmd.getOptionValue("fpType");
      boolean updateDictionaryFile = false;
      boolean hashUnknownFrag = false;
      Fingerprinter fprinter = Fingerprinter.createFingerprinter(type, 
                                          updateDictionaryFile, hashUnknownFrag);
      OEMolBase mol = new OEGraphMol();

      String inFile  = cmd.getOptionValue("i");
      oemolistream ifs = new oemolistream(inFile);
      
      double fract = 2D;
      String tmp = cmd.getOptionValue("sampleFract");
      if(tmp != null) fract = Double.parseDouble(tmp);
      
      Random rnd = new Random();
      
      LearningStrcutureCodeMapper mapper =
                           (LearningStrcutureCodeMapper) fprinter.getMapper(); 
      int dictSize = mapper.getMaxIdx()+1;
      int[] freq = new int[dictSize];
      
      while(oechem.OEReadMolecule(ifs, mol))
      {  iCounter++;
         if(rnd.nextDouble()< fract)
         {  fpCounter++;
            
            Fingerprint fp = fprinter.getFingerprint(mol);
            for(int bit : fp.getBits())
               freq[bit]++;
         }
         
         if(iCounter % 100 == 0) System.err.print(".");
         if(iCounter % 4000 == 0)
         {  System.err.printf( " %d %d %dsec\n",
                  iCounter, fpCounter, (System.currentTimeMillis()-start)/1000);
         }
      }
      
      System.err.printf("FPDictionarySorter: Read %d structures calculated %d fprints in %d sec\n",
            iCounter, fpCounter, (System.currentTimeMillis()-start)/1000);
      
      
      mapper.reSortDictionary(freq);
      mapper.writeDictionary();
      fprinter.close();
   }
   
   private static void exitWithHelp(Options options) {
      HelpFormatter formatter = new HelpFormatter();
      formatter.printHelp( "FPDictionarySorter", options );
      System.exit(1);
   }
}
