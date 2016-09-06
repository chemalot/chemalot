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

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.PosixParser;

import com.aestel.chemistry.openEye.SimComparatorFactory;
import com.aestel.chemistry.openEye.SphereExclusion;
import com.aestel.chemistry.openEye.fp.FPComparator;
import com.aestel.chemistry.openEye.fp.FPComparatorFact;

/**
 * SphereExclusion method implemented using internal fp toolkit and OEchem sdf reader.
 *
 * @author albertgo
 *
 */
public class SDFFPSphereExclusion
{
   private SDFFPSphereExclusion()
   {
   }


   public static void main(String...args) throws IOException
   {  // create command line Options object
      Options options = new Options();
      Option opt = new Option("in",true, "input file [.sdf,...]");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("out",true, "output file oe-supported");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("ref",true, "refrence file to be loaded before starting");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("fpTag",true, "field containing fingerpPrint");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("maxTanimoto",false, "If given the modified maxTanimoto will be used = common/(2*Max(na,nb)-common).");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("radius",true,
         "radius of exclusion sphere, exclude anything with similarity >= radius.");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("printSphereMatchCount",false, "check and print membership of candidates not "
            + " only to the first centroid which has sim >= radius but to all centroids"
            + " found up to that input. This will output a candidate multiple times."
            + " Implies checkSpheresInOrder.");
      options.addOption(opt);

      opt = new Option("checkSpheresInOrder",false, "For each candiate: compare to centroids from first to last (default is last to first)");
      options.addOption(opt);

      opt = new Option("printAll",false, "print all molecule, check includeIdx tag");
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

      // the only reason not to match centroids in reverse order id if
      // a non-centroid is to be assigned to multiple centroids
      boolean printSphereMatchCount = cmd.hasOption("printSphereMatchCount");
      boolean reverseMatch = ! cmd.hasOption("checkSpheresInOrder") && ! printSphereMatchCount;
      boolean printAll = cmd.hasOption("printAll") || printSphereMatchCount;
      boolean doMaxTanimoto = cmd.hasOption("maxTanimoto");
      String fpTag = cmd.getOptionValue("fpTag");
      double radius = Double.parseDouble(cmd.getOptionValue("radius"));
      String inFile  = cmd.getOptionValue("in");
      String outFile = cmd.getOptionValue("out");
      String refFile = cmd.getOptionValue("ref");

      SimComparatorFactory<OEMolBase, FPComparator, FPComparator> compFact
         = new FPComparatorFact(doMaxTanimoto, fpTag);
      SphereExclusion<FPComparator, FPComparator> alg =
         new SphereExclusion<FPComparator, FPComparator>(compFact, refFile,
                     outFile, radius, reverseMatch, printSphereMatchCount, printAll);
      alg.run(inFile);
      alg.close();
   }

   private static void exitWithHelp(Options options) {
      HelpFormatter formatter = new HelpFormatter();
      formatter.printHelp( "SphereExclusion", options );
      System.exit(1);
   }
}

