/*
   Copyright 2008-2014 Genentech Inc.

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
import java.io.IOException;
import java.io.InputStreamReader;

import openeye.oechem.OEExprOpts;
import openeye.oechem.OEMolBase;

import org.apache.commons.cli.*;

import com.aestel.chemistry.openEye.*;
import com.aestel.chemistry.openEye.nn.*;
import com.genentech.chemistry.openEye.*;
import com.genentech.chemistry.openEye.AAPathComparatorFact.AAPathCompareType;

/**
 * Command line program to compute a the nearest neighbors of a input file of compounds
 * using genetech similarity methods.
 *
 * This class only parses the command line options. the actual implementation is
 * performed by one of the implementations of {@link AbstractNNFinder} or
 * {@link AbstractNNMatrixFinder}.
 *
 * {@link AbstractNNMatrixFinder} is used if only one file is supplied and performance
 * an all by all comparison.
 *
 *  {@link AbstractNNFinder} is used if additionally a reference fiel is supplied
 * and computes a n*m comparison with n being the reference compoudns and m the
 * input compounds.
 *
 * @author albertgo
 *
 */
public class SDFMCSSNNFinder
{  static final Options options = new Options();
   static
   {   // create command line Options object

      Option opt = new Option("in",true, "input file [.sdf,...]");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("out",true, "output file oe-supported");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("ref",true, "refrence file to be loaded before starting, default compare to input");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("printAll",false, "Ouput each input record even when no neighbors are found.");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("idTag",true, "field containing id in reference file, used to create NNId fields. Only if maxNeighbors = 1 and minSim = 0");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("maxNeighbors",true, "maximum number of neighbors to be kept, default = 1");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("minSimilarity",true, "minimum similarity for molecule to be conisdered a neighbor, default = 0. Note for Cliff calculations this is the minimum cliff size.");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("countSimilarAbove",true, "Count number of similar at or above given similarity.");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("AAPathSim",true, "FAST|INTERMEDIATE|DEFAULT|FUZZY|FUZZYINTERMEDIATE|FUZZYFAST Use atom atom path match similarity of given type. Append version number (1-8) default=DEFAULT"
                                         + AAPathComparatorFact.DEFAULTVersion);
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("MCSSSimType",true, "DEFAULT|QueryRatio Compute MCSS and convert to sim using:. default=DEFAULT\n"
                      + "   2 * Default = nAtMatch/(nAtQuery+nAtCand - nAtMatch) + nBdMatch/(nBdQuery+nBdCand - nBdMatch)\n" 
                      + "   QueryRatio  = nAtMatch/nAtQuery");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("atomMatch", true, "For MCSS only: Sequence of none|default|hcount|noAromatic specifing how atoms are matched cf. oe document.\n"
            +"noAromatic can be used to make terminal atoms match aliphatic and aromatic atoms.\n"
            +"Queryfeatures are considered only if default is used.");
      options.addOption(opt);
      
      opt = new Option("bondMatch", true, "For MCSS only: Sequence of none|default specifing how bonds are matched cf. oe document.");
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

      opt = new Option("cliffPropertyTag",true, "If this is given then the activity cliff is computed instead of the similarity.");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("cliffSimAt2Delta",true, "Use alternate cliff formula:\n"
                                               +"cliff = abs(delta(prop)) * e^(a*sim^b)\n"
                                               +"a = ln(foldAt1Sim)\n"
                                               +"b = log(ln(2)/a), simAt2Delta)\n"
                                               +"simAt2Delta: similarity at which delta(prop) is doubled.\n"
                                               +"foldAt1Sim:  Maximum Cliff size (sim~=1) is foldAt1Sim * delta(props)");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("cliffFoldAt1Sim",true, "cf. -cliffSimAt2Delta");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("cliffPropertyLogarithmic",false, "If this is given then the lgo10 of the property is taken when computing the cliff value.");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("cliffIgnoreModifier",false, "Ignore '<>= ~' infront of activity value.");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("cliffMinSimilarity",true, "Minimum similarity between cliff pairs (default 0)");
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

      String idTag = cmd.getOptionValue("idTag");
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
         exitWithHelp("-outputDuplicates will not work with outputVTab");
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
         System.err.println("WARNING: printAll ignored for tab output!\n");

      SimComparatorFactory<OEMolBase, OEMolBase, SimComparator<OEMolBase>> compFact;
      compFact = getComparatorFactory(cmd);

      if( refFile == null )
      {  // no reference file; run all by all comparison
         performMatrixNNSearch(inFile, outFile, tabOutput, compFact, minSim,
                  maxNeighbors, idTag, nCpu, countAboveSimilarityStr, printAll);

      }else
      {  // refrence file; compare inFile to refFile
         performReferenceSearch(inFile, refFile, outFile, tabOutput, compFact,
               minSim, maxNeighbors, idTag, nCpu, countAboveSimilarityStr,
               outputDuplicates, printAll);
      }

   }


   /** parse options and get similarity comparator factory */
   private static SimComparatorFactory<OEMolBase, OEMolBase, SimComparator<OEMolBase>>
                     getComparatorFactory(CommandLine cmd)
   {  SimComparatorFactory<OEMolBase, OEMolBase, SimComparator<OEMolBase>> compFact;
      if( cmd.hasOption("AAPathSim") )
      {  if( cmd.hasOption("MCSSSimType") )
            throw new Error("-AAPathSim and -MCSSSimType may not be used together");
      
         if( cmd.hasOption("atomMatch") || cmd.hasOption("bondMatch") )
            throw new Error("-AAPathSim does not support '-atomMatch' or '-bondMatch'");
         
         String aaPathSimType = cmd.getOptionValue("AAPathSim");
         int version = AAPathComparatorFact.DEFAULTVersion;  // current default version
         if( aaPathSimType.matches(".*\\d") )
         {  version = aaPathSimType.charAt(aaPathSimType.length()-1) - '0';
            aaPathSimType = aaPathSimType.substring(0,aaPathSimType.length()-1);
         }
         AAPathCompareType type = AAPathCompareType.valueOf(aaPathSimType);
         compFact = new AAPathComparatorFact(type, version);
      }
      else
      {  int atExpr = OEExprOpts.DefaultAtoms;
         String atomMatch = cmd.getOptionValue("atomMatch");
         if( atomMatch == null ) atomMatch = "";
         atomMatch = '|' + atomMatch.toLowerCase() + '|';
         if( atomMatch.startsWith("|none") )      atExpr = 0;
         if( atomMatch.contains("|hcount|") )     atExpr |= OEExprOpts.HCount;
         if( atomMatch.contains("|noAromatic|") ) atExpr &= (~ OEExprOpts.Aromaticity);
   
         int bdExpr = OEExprOpts.DefaultBonds;
         String bondMatch = cmd.getOptionValue("bondMatch");
         if( bondMatch == null ) bondMatch = "";
         bondMatch = '|' + bondMatch.toLowerCase() + '|';
         if( bondMatch.startsWith("|none") ) bdExpr = 0;
         
         MCSSCompareType mcssType = MCSSCompareType.DEFAULT;
         if( cmd.hasOption("MCSSSimType") )
            mcssType = MCSSCompareType.toEnum(cmd.getOptionValue("MCSSSimType"));

         compFact = new MCSSComparatorFact(mcssType, atExpr, bdExpr);
      }

      if( cmd.hasOption("cliffPropertyTag") )
      {  compFact = getCliffComparatorFact(compFact, cmd);
      }
      return compFact;
   }


   /** construct a CliffComaprator as requested by the command line options */
   private static SimComparatorFactory<OEMolBase, OEMolBase, SimComparator<OEMolBase>>
                     getCliffComparatorFact( SimComparatorFactory<OEMolBase, OEMolBase,
                                             SimComparator<OEMolBase>> compFact,
                                             CommandLine cmd)
   {  String d;
      double cliffFoldAt1Sim = 0D;
      double cliffSimAt2Delta = 0D;

      d = cmd.getOptionValue("cliffSimAt2Delta");
      if( d != null ) cliffSimAt2Delta = Double.parseDouble(d);

      d = cmd.getOptionValue("cliffFoldAt1Sim");
      if( d != null ) cliffFoldAt1Sim = Double.parseDouble(d);

      if( (cliffSimAt2Delta <= 0D || cliffSimAt2Delta >= 1D || cliffFoldAt1Sim <= 2D)
               && ! (cliffSimAt2Delta == 0 && cliffFoldAt1Sim == 0))
         exitWithHelp("-cliffSimAt2Delta must be between 0 and 1 (suggested .4)\n"
                     +"-cliffFoldAt1Sim must be > 2 (suggested 100)");

      double minCliffSim = 0D;
      String cliffPropTag = cmd.getOptionValue("cliffPropertyTag");

      boolean cliffPropLog = cmd.hasOption("cliffPropertyLogarithmic");

      boolean cliffIgnoreModifier = cmd.hasOption("cliffIgnoreModifier");

      d = cmd.getOptionValue("cliffMinSimilarity");
      if( d != null ) minCliffSim = Double.parseDouble(d);

      if( cliffPropLog && cliffPropTag == null )
         exitWithHelp("cliffPropertyLogarithmic needs cliffPropertyTag option");

      if( cliffIgnoreModifier && cliffPropTag == null )
         exitWithHelp("cliffIgnoreModifier needs cliffPropertyTag option");

      return new CliffComparatorFact(compFact, cliffPropTag, cliffIgnoreModifier, cliffPropLog,
                                     minCliffSim, cliffFoldAt1Sim, cliffSimAt2Delta);
   }


   /** compare compounds in inFile to eachother */
   private static void performMatrixNNSearch(
            String inFile, String outFile, String tabOutput,
            SimComparatorFactory<OEMolBase, OEMolBase, SimComparator<OEMolBase>> compFact,
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
            c = new MultiNNMatrixFinderVTConsumer(outFile);
         else if( "tab".equalsIgnoreCase(tabOutput) )
            c = new MultiNNMatrixFinderTabConsumer(outFile, idTag);
         else
            c = new MultiNNMatrixFinderConsumer(outFile, countAboveSimilarityStr);

         alg = new MultiNNMatrixFinder<OEMolBase, SimComparator<OEMolBase>>(
                        inFile, c, compFact, maxNeighbors, minSim, printAll, countAboveSimilarity);

         if( "tab".equalsIgnoreCase(tabOutput) )
            ((MultiNNMatrixFinderTabConsumer)c).setMatrixSize(alg.getObjectCount());
      }else
      {  NNMatrixFinderConsumerInterface c;
         if( "vTab".equalsIgnoreCase(tabOutput) )
            c = new NNMatrixFinderVTConsumer(outFile);
         else
            c = new NNMatrixFinderConsumer(outFile, countAboveSimilarityStr);

         alg = new NNMatrixFinder<OEMolBase, SimComparator<OEMolBase>>(
                                                   inFile, c, compFact, countAboveSimilarity);
      }
      MultiThreadMatrixRunner<OEMolBase, MCSSComparator> runner
         = new MultiThreadMatrixRunner<OEMolBase, MCSSComparator>(alg, nCpu);
      runner.run();
      runner.close();
   }


   /**
    * For each compound in inFile comapre to each in refFile.
    */
   private static void performReferenceSearch(
            String inFile, String refFile, String outFile, String tabOutput,
            SimComparatorFactory<OEMolBase, OEMolBase, SimComparator<OEMolBase>> compFact,
            double minSim, int maxNeighbors, String idTag, int nCpu,
            String countAboveSimilarityStr, boolean outputDuplicates, boolean printAll)
   throws IOException
   {  double countAboveSimilarity = Double.MAX_VALUE;
      if( countAboveSimilarityStr != null )
         countAboveSimilarity = Double.parseDouble(countAboveSimilarityStr);

      MultiThreadAlgortihm nnAlg;

      if( maxNeighbors > 1 || minSim > 0 )  // use consumer to output n NN's
      {  NNMultiFinderConsumerInterface c;
         if( "vTab".equalsIgnoreCase(tabOutput) )
            c = new NNMultiFinderVTConsumer(outFile, idTag != null);
         else if( outputDuplicates )
            c = new NNMultiFinderDuplConsumer(outFile, idTag != null, countAboveSimilarityStr);
         else
            c = new NNMultiFinderConsumer(outFile, idTag != null, countAboveSimilarityStr);

         nnAlg = new MultiNNFinder<OEMolBase, SimComparator<OEMolBase>>(inFile, c,
                  compFact, refFile, idTag, maxNeighbors, minSim, printAll, countAboveSimilarity);

      } else
      {  // use consumer to output only single NN
         NNFinderConsumerInterface c;
         if( "vTab".equalsIgnoreCase(tabOutput) )
            c = new NNFinderVTConsumer(outFile, idTag);
         else
            c = new NNFinderConsumer(outFile, idTag, countAboveSimilarityStr);

         nnAlg = new NNFinder<OEMolBase, SimComparator<OEMolBase>>(
                                          inFile, c, compFact,refFile, idTag, countAboveSimilarity );
      }

      MultiThreadRunner runner = new MultiThreadRunner(nnAlg, nCpu);
//System.err.print("Waiting, hit enter: ");System.in.read();
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
         +"       is given.\n"
         +"    If no refFile is given an all to all comparison is done.\n"
         +"  The following fields will be added to your output file:\n"
         +"    idx the index in the input file (only if no refFile given, because\n"
         +"        the order of the compounds is not preserved\n"
         +"    NNIdx index ',' separated of near neighbors\n"
         +"    NNSim similarity ',' separated of near neighbors\n";
      formatter.printHelp( 120, "sdfMCSSNNFinder", explain, options, "" );
      System.exit(1);
   }


   private SDFMCSSNNFinder()
   {  // singleton
   }
}

