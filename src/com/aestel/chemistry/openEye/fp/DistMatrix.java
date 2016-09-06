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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.PosixParser;

/**
 * Read tab separated file with "id\tFingerPrint" and output
 * Similarity matrix for the comparison of each fingerprint with each other.
 * Output result as tab separated file.
 * @author albertgo
 *
 */
public class DistMatrix
{
   public static void main(String...args) throws IOException
   {  long start = System.currentTimeMillis();

      // create command line Options object
      Options options = new Options();
      Option opt = new Option("i",true, "input file [.tsv from FingerPrinter]");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("o",true, "outpur file [.tsv ");
      opt.setRequired(true);
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

      if( args.length != 0 ) exitWithHelp(options);

      String file = cmd.getOptionValue("i");
      BufferedReader in = new BufferedReader(new FileReader(file));

      file = cmd.getOptionValue("o");
      PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(file)));

      ArrayList<Fingerprint> fps = new ArrayList<Fingerprint>();
      ArrayList<String> ids = new ArrayList<String>();
      String line;
      while( (line = in.readLine()) != null )
      {  String[] parts = line.split("\t");
         if(parts.length == 3)
         {  ids.add(parts[0]);
            fps.add(new ByteFingerprint(parts[2]));
         }
      }
      in.close();

      out.print("ID");
      for(int i=0; i< ids.size(); i++)
      {  out.print('\t');
         out.print(ids.get(i));
      }
      out.println();

      for(int i=0; i<ids.size(); i++)
      {  out.print(ids.get(i));
         Fingerprint fp1 = fps.get(i);

         for(int j=0; j<=i;j++)
         {  out.printf("\t%.4g", fp1.tanimoto( fps.get(j)));
         }
         out.println();
      }
      out.close();

      System.err.printf("Done %d fingerprints in %.2gsec\n",
            fps.size(), (System.currentTimeMillis() - start)/1000D);
   }

   private static void exitWithHelp(Options options)
   {  HelpFormatter formatter = new HelpFormatter();
      formatter.printHelp("DistMatrix", options);
      System.exit(1);
   }

}
