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
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.regex.Pattern;

import openeye.oechem.OEGraphMol;
import openeye.oechem.OEMolBase;
import openeye.oechem.oechem;
import openeye.oechem.oemolithread;
import openeye.oechem.oemolothread;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.PosixParser;

import com.genentech.chemistry.openEye.MDLSSSMatcher;
/**
 * @param <T> comparable
 * @param <Y> comparator
 */
public class SDFMDLSSSMatcher
{  private final MDLSSSMatcher matcher;
   private final boolean firstMatch;
   private final ExecutorService executor;
   private final ExecutorCompletionService<Boolean> completionService;
   private final oemolothread ofs;
   private final int nCpu;
   private oemolithread ifs;
   private MDLSSSFind[] finder;
   public boolean printAll;

   public SDFMDLSSSMatcher(String qFile, String outFile, boolean firstMatch, boolean printAll, int nCpu)
   {  this.firstMatch = firstMatch;
      this.printAll   = printAll;
      this.matcher = new MDLSSSMatcher(qFile);
      this.ofs     = new oemolothread(outFile);
      this.nCpu    = nCpu;

      this.executor = Executors.newFixedThreadPool(nCpu);
      this.completionService = new ExecutorCompletionService<Boolean>(executor);

      finder = new MDLSSSFind[nCpu];
   }

   void run(String inFile) throws InterruptedException
   {  long start = System.currentTimeMillis();
      ifs = new oemolithread(inFile);

      // start nCpu instances of MDLSSSFind and submit them to the completionService
      for(int i =0; i< nCpu; i++)
      {  finder[i] = new MDLSSSFind();
         completionService.submit(finder[i]);
      }

      int iCount = 0;
      int lastICount = 0;
      int nRunning = nCpu;
      // poll the running threads and output progress until all are done
      while( nRunning > 0 )
      {  if( completionService.poll(2, TimeUnit.SECONDS) != null)
         {  nRunning--;
         }

         // print status report to stderr
         iCount = 0;
         for( int i=0; i<nCpu; i++)
            iCount += finder[i].count;

         if( iCount == lastICount ) continue;

         if( iCount / (4000*nCpu) > lastICount / (4000*nCpu))
         {  System.err.printf( " %d %dsec\n",
                  iCount, (System.currentTimeMillis()-start)/1000);
         }else
         {  for( int i=lastICount; i< iCount - (100*nCpu); i+=(100*nCpu))
               System.err.print(".");
         }
         lastICount = iCount;
      }

      ifs.close();
      ifs.delete();
      ifs = null;

      int nMatches = 0;
      for( int i=0; i<nCpu; i++)
         nMatches += finder[i].nMatches;

      inFile = inFile.replaceAll(".*" + Pattern.quote(File.separator), "");
      System.err.printf("\nSDFMDLSSSMatcher: Read %d structures from %s, found %d matches. nCpu=%d %d sec\n",
            iCount, inFile, nMatches, nCpu, (System.currentTimeMillis()-start)/1000);
   }

   class MDLSSSFind implements Callable<Boolean>
   {  volatile int count = 0;
      int nMatches = 0;

      @Override
      public Boolean call()
      {  OEMolBase mol = new OEGraphMol();
         while( oechem.OEReadMolecule(ifs, mol) )
         {  if( matcher.findMatches(mol, firstMatch))
            {  oechem.OEWriteMolecule(ofs, mol);
               nMatches++;
            }else if( printAll )
            {  oechem.OEWriteMolecule(ofs, mol);
            }

            count++;
         }
         mol.delete();

         return Boolean.FALSE;
      }
   }

   public void close()
   {  matcher.close();
      executor.shutdown();
      ofs.close();
   }



   public static void main(String...args) throws IOException, InterruptedException
   {  oechem.OEUseJavaHeap(false);

      // create command line Options object
      Options options = new Options();
      Option opt = new Option("in",true, "input file [.sdf,...]");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("out",true, "output file oe-supported");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("ref",true, "refrence file with MDL query molecules");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("anyMatch",false, "if set all matches are reported not just the first.");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("printAll",false, "if set even compounds that do not macht are outputted.");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("nCpu",true, "number of CPU's used in parallel, dafault 1");
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
      if( args.length > 0 )
      {  exitWithHelp("Unknown param: " + args[0], options);
      }

      if(cmd.hasOption("d"))
      {  System.err.println("Start debugger and press return:");
         new BufferedReader(new InputStreamReader(System.in)).readLine();
      }

      int nCpu     = 1;
      boolean firstMatch = ! cmd.hasOption("anyMatch");
      boolean printAll = cmd.hasOption("printAll");


      String d     = cmd.getOptionValue("nCpu");
      if( d != null ) nCpu = Integer.parseInt(d);

      String inFile  = cmd.getOptionValue("in");
      String outFile = cmd.getOptionValue("out");
      String refFile = cmd.getOptionValue("ref");

      SDFMDLSSSMatcher matcher = new SDFMDLSSSMatcher(refFile, outFile, firstMatch, printAll, nCpu);
      matcher.run(inFile);
      matcher.close();
   }

   private static void exitWithHelp(String msg , Options options)
   {  System.err.println( msg );
      exitWithHelp(options);
   }


   private static void exitWithHelp(Options options) {
      HelpFormatter formatter = new HelpFormatter();
      String explain = "Will match each of the queries in the reference file against\n"
         +"  all molecule in the input file. The queries in the reference file are\n"
         +"  interpreted as MDLQuery molecules.\n"
         +"  Only compound matching at least one reference query will be outputted.\n"
         +"  These molecules will contain the additional SSSMatchName tag with the\n"
         +"  value of the molfile-title of the matching query molecules.\n";
      formatter.printHelp( 120, "SDFSSSMatcher", explain, options, "" );
      System.exit(1);
   }


}
