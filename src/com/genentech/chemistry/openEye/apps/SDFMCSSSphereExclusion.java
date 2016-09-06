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

import openeye.oechem.OEMolBase;

import org.apache.commons.cli.*;

import com.aestel.chemistry.openEye.SimComparator;
import com.aestel.chemistry.openEye.SimComparatorFactory;
import com.aestel.chemistry.openEye.SphereExclusion;
import com.genentech.chemistry.openEye.AAPathComparatorFact;
import com.genentech.chemistry.openEye.MCSSComparatorFact;
import com.genentech.chemistry.openEye.MCSSCompareType;
import com.genentech.chemistry.openEye.AAPathComparatorFact.AAPathCompareType;

/**
 * Command line program to diversity selection usign the SphereExclusion algorithm
 * with genetech similarity methods.
 *
 * This class only parses the command line options. the actual implementation is
 * performed in {@link SphereExclusion}
 *
 * @author albertgo
 *
 */
public class SDFMCSSSphereExclusion
{  private static final Options options = new Options();
   static
   {  // create command line Options object
      Option opt = new Option("in",true, "candidate file: input file oe-supported");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("out",true, "output file oe-supported");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("ref",true, "refrence file to be loaded before starting");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("radius",true,
         "radius of exclusion sphere, exclude anything with similarity >= radius.");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("AAPathSim",true, "FAST|INTERMEDIATE|DEFAULT|FUZZY|FUZZYINTERMEDIATE|FUZZYFAST Use atom atom path match similarity of given type. Append version number (1-8) default=DEFAULT"
               + AAPathComparatorFact.DEFAULTVersion);
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("printSphereMatchCount",false, "Assign non-centroids to all spheres with sim>= radius");
      options.addOption(opt);

      opt = new Option("checkSpheresInOrder",false, "For each candiate: compare to centroids from first to last (default is last to first)");
      options.addOption(opt);

      opt = new Option("printAll",false, "print all molecule, check includeIdx tag");
      options.addOption(opt);
   }


   public static void main(String...args) throws IOException
   {  CommandLineParser parser = new PosixParser();
      CommandLine cmd = null;
      try
      {  cmd = parser.parse( options, args);
      } catch(Exception e)
      {  System.err.println(e.getMessage());
         exitWithHelp();
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
      double radius = Double.parseDouble(cmd.getOptionValue("radius"));
      String inFile  = cmd.getOptionValue("in");
      String outFile = cmd.getOptionValue("out");
      String refFile = cmd.getOptionValue("ref");

      SimComparatorFactory<OEMolBase, OEMolBase, SimComparator<OEMolBase>> compFact;
      compFact = getComparatorFactory(cmd);

      SphereExclusion<OEMolBase, SimComparator<OEMolBase>> alg =
         new SphereExclusion<OEMolBase, SimComparator<OEMolBase>>(compFact, refFile,
                  outFile, radius, reverseMatch, printSphereMatchCount, printAll);
      alg.run(inFile);
      alg.close();
   }

   /**
    * Construct the ComparatorFactory from the command line options
    */
   private static SimComparatorFactory<OEMolBase, OEMolBase, SimComparator<OEMolBase>>
                     getComparatorFactory(CommandLine cmd)
   {
      SimComparatorFactory<OEMolBase, OEMolBase, SimComparator<OEMolBase>> compFact;
      if( cmd.hasOption("AAPathSim") )
      {  String aaPathSimType = cmd.getOptionValue("AAPathSim");
         int version = AAPathComparatorFact.DEFAULTVersion;
         if( aaPathSimType.matches(".*\\d") )
         {  version = aaPathSimType.charAt(aaPathSimType.length()-1) - '0';
            aaPathSimType = aaPathSimType.substring(0,aaPathSimType.length()-1);
         }
         AAPathCompareType type = AAPathCompareType.valueOf(aaPathSimType);
         compFact = new AAPathComparatorFact(type, version);
      }
      else
      {  compFact = new MCSSComparatorFact(MCSSCompareType.DEFAULT);
      }
      return compFact;
   }

   private static void exitWithHelp()
   {  HelpFormatter formatter = new HelpFormatter();
      formatter.printHelp( "SDFMCSSSphereExclusion", options );
      System.exit(1);
   }


   private SDFMCSSSphereExclusion()
   {  // singlton
   }
}

