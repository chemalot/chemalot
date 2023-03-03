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

package com.aestel.chemistry.openEye.fp.apps;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

import openeye.oechem.OEMolBase;

import org.apache.commons.cli.*;

import com.aestel.chemistry.openEye.*;
import com.aestel.chemistry.openEye.fp.FPComparator;
import com.aestel.chemistry.openEye.fp.FPComparatorFact;
import com.aestel.chemistry.openEye.nn.*;

/**
 * SphereExclusion method implemented using internal fp toolkit and OEchem sdf reader.
 *
 * @author albertgo
 *
 */
public class SDFFPNNFinder
{  private static final Options options = new Options();
   static
   {  // create command line Options object
      Option opt = new Option("in",true, "input file [.sdf,...]");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("out",true, "output file oe-supported");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("ref",true, "refrence file to be loaded before starting, default compare to input");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("maxNeighbors",true, "maximum number of neighbors to be kept, default = 1");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("minSimilarity",true, "minimum similarity for molecule to be conisdered a neighbor, default = 0");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("countSimilarAbove",true, "Count number of similar at or above given similarity.");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("printAll",false, "Ouput each input record even when no neighbors are found.");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("fpTag",true, "field containing fingerpPrint");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("idTag",true, "field containing id in reference file, used to create NNId fields. Only if maxNeighbors = 1 and minSim = 0");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("maxTanimoto",false, "If given the modified maxTanimoto will be used = common/(2*Max(na,nb)-common).");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("outputDuplicates",false, "If a candidate matches multiple reference compounds output the candidte once for each match.");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("tabOutput", true, "tab|vTab output as table NxN or vertical table with smiles.smiles<tab>NNSim");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("nCpu",true, "number of CPU's used in parallel, dafault 1");
      opt.setRequired(false);
      options.addOption(opt);
   }



   public static void main(String...args) throws IOException
   {
      CommandLineParser parser = new PosixParser();
      CommandLine cmd = null;
      try
      {  cmd = parser.parse( options, args);
      } catch(Exception e)
      {  System.err.println(e.getMessage());
         exitWithHelp();
      }

      args = cmd.getArgs();
      if( args.length > 0 )
      {  exitWithHelp("Unknown param: " + args[0]);
      }

      if(cmd.hasOption("d"))
      {  System.err.println("Start debugger and press return:");
         new BufferedReader(new InputStreamReader(System.in)).readLine();
      }

      int nCpu     = 1;
      int maxNeighbors = 1;
      double minSim = 0D;

      String fpTag = cmd.getOptionValue("fpTag");
      String idTag = cmd.getOptionValue("idTag");
      boolean doMaxTanimoto = cmd.hasOption("maxTanimoto");
      boolean printAll = cmd.hasOption("printAll");

      String d     = cmd.getOptionValue("nCpu");
      if( d != null ) nCpu = Integer.parseInt(d);

      d            = cmd.getOptionValue("maxNeighbors");
      if( d != null ) maxNeighbors = Integer.parseInt(d);

      d            = cmd.getOptionValue("minSimilarity");
      if( d != null ) minSim = Double.parseDouble(d);

      String countAboveSimilarityStr = cmd.getOptionValue("countSimilarAbove");

      String inFile  = cmd.getOptionValue("in");
      String outFile = cmd.getOptionValue("out");
      String refFile = cmd.getOptionValue("ref");

      String tabOutput = cmd.getOptionValue("tabOutput");
      boolean outputDuplicates = cmd.hasOption("outputDuplicates");

      if( outputDuplicates && tabOutput != null )
         exitWithHelp("-outputDuplicates will not work with tabOutput");

      if( outputDuplicates && refFile == null )
         exitWithHelp("-outputDuplicates requires -ref ");

      if( "tab".equalsIgnoreCase(tabOutput) && refFile != null )
         exitWithHelp("-tabOutput tab: does not work with reference file");

      if( "tab".equalsIgnoreCase(tabOutput) && maxNeighbors == 1 )
         exitWithHelp("-tabOutput tab: does not make sense with -maxNeighbors = 1");

      if( cmd.hasOption("countSimilarAbove") && tabOutput != null )
         exitWithHelp("-countSimilarAbove not supported for tab or vTab output");

      if( printAll && ! (maxNeighbors > 1 || minSim > 0) )
         exitWithHelp("printAll only supported if: maxNeighbors > 1 or minSim > 0");

      if( printAll && tabOutput != null )
         System.err.println("WARNING: printAll ignored tor tab output!\n");


      SimComparatorFactory<OEMolBase, FPComparator, FPComparator> compFact
         = new FPComparatorFact(doMaxTanimoto, fpTag);

      if( refFile == null )
      {  perfromMatrixNNSearch(inFile, outFile, tabOutput, compFact, minSim,
               maxNeighbors, idTag, nCpu, countAboveSimilarityStr, printAll);

      }else
      {  performReferenceSearch(inFile, refFile, outFile, tabOutput, compFact,
               minSim, maxNeighbors, idTag, nCpu, countAboveSimilarityStr,
               outputDuplicates, printAll);
      }

   }


   /** compare compounds in inFile to eachother */
   private static void perfromMatrixNNSearch(
            String inFile, String outFile, String tabOutput,
            SimComparatorFactory<OEMolBase, FPComparator, FPComparator> compFact,
            double minSim, int maxNeighbors, String idTag, int nCpu,
            String countAboveSimilarityStr, boolean printAll)
   throws IOException
   {  double countAboveSimilarity = Double.MAX_VALUE;
      if( countAboveSimilarityStr != null )
         countAboveSimilarity = Double.parseDouble(countAboveSimilarityStr);

      MultiThreadMatrixAlgortihm alg;

      if( maxNeighbors > 1 || minSim > 0 )
      {  MultiNNMatrixFinderConsumerInterface c;
         if( "vTab".equalsIgnoreCase(tabOutput) )
            c = new MultiNNMatrixFinderVTConsumer(outFile, idTag);
         else if( "tab".equalsIgnoreCase(tabOutput) )
            c = new MultiNNMatrixFinderTabConsumer(outFile, idTag);
         else
            c = new MultiNNMatrixFinderConsumer(outFile, countAboveSimilarityStr);

         alg = new MultiNNMatrixFinder<FPComparator, FPComparator>(inFile, c,
                           compFact, maxNeighbors, minSim, printAll, countAboveSimilarity);

         if( "tab".equalsIgnoreCase(tabOutput) )
            ((MultiNNMatrixFinderTabConsumer)c).setMatrixSize(alg.getObjectCount());
      }else
      {  NNMatrixFinderConsumerInterface c;
         if( "vTab".equalsIgnoreCase(tabOutput) )
            c = new NNMatrixFinderVTConsumer(outFile, idTag);
         else
            c = new NNMatrixFinderConsumer(outFile, countAboveSimilarityStr);

         alg = new NNMatrixFinder<FPComparator, FPComparator>(
                                          inFile, c, compFact, countAboveSimilarity);
      }

      MultiThreadMatrixRunner<FPComparator, FPComparator> runner
         = new MultiThreadMatrixRunner<FPComparator, FPComparator>(alg, nCpu);
      runner.run();
      runner.close();
   }


   /**
    * For each compound in inFile comapre to each in refFile.
    */
   private static void performReferenceSearch(
            String inFile, String refFile, String outFile, String tabOutput,
            SimComparatorFactory<OEMolBase, FPComparator, FPComparator> compFact,
            double minSim, int maxNeighbors, String idTag, int nCpu,
            String countAboveSimilarityStr, boolean outputDuplicates, boolean printAll)
   throws IOException
   {  double countAboveSimilarity = Double.MAX_VALUE;
      if( countAboveSimilarityStr != null )
         countAboveSimilarity = Double.parseDouble(countAboveSimilarityStr);

      MultiThreadAlgortihm nnAlg;

      if( maxNeighbors > 1 || minSim > 0 )
      {  NNMultiFinderConsumerInterface c;
         if( "vTab".equalsIgnoreCase(tabOutput) )
            c = new NNMultiFinderVTConsumer(outFile, idTag, idTag != null);
         else if( outputDuplicates )
            c = new NNMultiFinderDuplConsumer(outFile, idTag != null, countAboveSimilarityStr);
         else
            c = new NNMultiFinderConsumer(outFile, idTag != null, countAboveSimilarityStr);

         nnAlg = new MultiNNFinder<FPComparator, FPComparator>( inFile, c,
                  compFact,refFile, idTag, maxNeighbors, minSim, printAll, countAboveSimilarity);

      } else
      {  NNFinderConsumerInterface c;
         if( "vTab".equalsIgnoreCase(tabOutput) )
            c = new NNFinderVTConsumer(outFile, idTag);
         else
            c = new NNFinderConsumer(outFile, idTag, countAboveSimilarityStr);

         nnAlg = new NNFinder<FPComparator, FPComparator>(
                                          inFile, c, compFact,refFile, idTag, countAboveSimilarity );
      }

      MultiThreadRunner runner = new MultiThreadRunner(nnAlg, nCpu);
      runner.run();
      runner.close();
   }


   private static void exitWithHelp(String msg)
   {  System.err.println( msg );
      exitWithHelp();
   }


   private static void exitWithHelp() {
      HelpFormatter formatter = new HelpFormatter();
      String explain = "Will calculate the nearest neighbors for the input compounds\n"
         +"  Two modes are supported:\n"
         +"    If a reference file is given the NN from to the compounds in the reference file\n"
         +"       is given. First the reference compounds are written to the output then the in-compounds.\n"
         +"    If no refFile is given an all to all comparison is done.\n"
         +"  The following fields will be added to your output file:\n"
         +"    idx the index in the input file (only if no refFile given, because\n"
         +"        the order of the compounds is not preserved\n"
         +"    NNIdx index ',' separated of near neighbors\n"
         +"    NNSim similarity ',' separated of near neighbors\n";
      formatter.printHelp( 120, "SDFNNFinder", explain, options, "" );
      System.exit(1);
   }


   private SDFFPNNFinder()
   {  // singlton
   }
}

