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
package com.genentech.chemistry.tool.sdfAggregator;

/**
 * This class provides the framework to call aggregation functions such as count, mean, etc.
 * @author Johnny Wu Aug 11, 2011
 *
 */

import java.io.IOException;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import openeye.oechem.*;

import org.apache.commons.cli.*;


public class SDFAggregator
{
   private static final String PACKAGEName;
   private static CommandLine cl;
   private static Options options;
   private static PosixParser parser;
   private static String inFile;
   private static String outFile;
   private static OUTPUTMODE outputMode;
   private static List<String> groups;
   private static AggInterface[] aAgg;
   private static String currentGroupValue;
   private static oemolostream ofs;
   private static oemolistream ifs;
   private static oemolostream oss;
   private static String[] functionArr = { "concatenate", "count", "max", "mean", 
         "min", "rank", "rankPct", "sum", "stddev", "stderr", "median", 
         "boltzmannAvg", "boltzmannProbability"};

   enum OUTPUTMODE
   {  FIRST, LAST, ALL }


   static
   {  String classname = SDFAggregator.class.getName();
      PACKAGEName = classname.substring(0, classname.lastIndexOf('.'));
   }

   /**
    * Generates help string of all aggregator functions
    *
    */
   private static String getAggregationFunctionDescriptions()
   {
      StringBuilder sb = new StringBuilder(2000);
      String desc = "";

      for (String funcName : functionArr)
      {  try
         {  desc = (String)Class.forName(PACKAGEName + "." + funcName)
                                    .getField("DESCRIPTION").get(null);
         } catch (NoSuchFieldException e)
         {
            System.err.println(e.getMessage());
            System.err.printf("Aggregation function %s does not implement DESCRIPION", funcName);
            exitWithHelp(options);
         } catch (IllegalAccessException e)
         {
            System.err.println(e.getMessage());
            System.err.println("Function does not exist: " + funcName);
            exitWithHelp(options);
         } catch (ClassNotFoundException e)
         {
            e.printStackTrace();
            System.err.println(e.getMessage());
            System.err.println("Function does not exist: " + funcName);
            exitWithHelp(options);
         }

         sb.append(desc);
      }
      return sb.toString();
   }

   private static void processOptions(String[] args)
   {  Option opt;
      options = new Options();
      parser = new PosixParser();

      opt = new Option("in", true,
               "Input file oe-supported Use .sdf to specify the file type.");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("out", true,
               "Output file oe-supported Use .sdf to specify the file type.");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option(
               "groupby",
               true,
               "Tags to group by (i.e. ClusterID, AssayName). Assumes tags are presorted (i.e. through SDFSorter). "
                        + "Use multiple -groupby tags to group by more than one field.");
      opt.setRequired(false);
      options.addOption(opt);
      groups = new ArrayList<String>();

      String funcListStr = getAggregationFunctionDescriptions();
      opt = new Option("function", true, "Aggregator function:\n" + funcListStr);
      opt.setRequired(true);
      options.addOption(opt);

      options.addOption("outputmode", true,
               "[all|first(def)|last] output function results to all or only first/last entry.");

      /*
       * private static void processOptions (String[] args){
       */
      try
      {
         cl = parser.parse(options, args);
      } catch (Exception exp)
      {  // catch (ParseException exp) {
         System.err.println(exp.getMessage());
         exitWithHelp(options);
      }

      // TODO: Check files exist
      inFile = cl.getOptionValue("in");
      outFile = cl.getOptionValue("out");

      outputMode = OUTPUTMODE.FIRST;
      if (cl.hasOption("outputmode"))
      {  String outModeStr = cl.getOptionValue("outputmode");
         if (outModeStr.matches("(?i)all"))
            outputMode = OUTPUTMODE.ALL;
         else if( outModeStr.equalsIgnoreCase("last"))
            outputMode = OUTPUTMODE.LAST;
      }

      String[] grpByArr = cl.getOptionValues("groupby");
      if( grpByArr != null )
         groups.addAll(Arrays.asList(grpByArr));

      String[] funcstr = cl.getOptionValues("function");

      aAgg = createAggregationFunction(funcstr);
   }


   private static AggInterface[] createAggregationFunction(String[] funcStrs)
   {  String outTag = null;
      String resFormat = null;
      String functionName = null;
      String fnargs = null;
      AggInterface[] ret = new AggInterface[funcStrs.length];

      for( int i=0; i<funcStrs.length; i++ )
      {  String funcstr = funcStrs[i];

         Pattern pat = Pattern.compile(" *([^:]+)(:.*)? *= *(\\w+) *\\( *(.*)\\) *");
         Matcher mat = pat.matcher(funcstr);
         if( mat.matches() )
         {  outTag = mat.group(1);
            resFormat = mat.group(2);
            functionName = mat.group(3);
            fnargs = mat.group(4);

            if( resFormat != null ) resFormat = resFormat.substring(1);

         } else // no output name given => construct from function name
         {  pat = Pattern.compile(" *(\\w+) *\\( *(.+) *\\) *");
            mat = pat.matcher(funcstr);

            if( mat.matches() )
            {  functionName = mat.group(2);
               fnargs = mat.group(3);
               outTag = functionName + '-' + fnargs;

            } else if( funcstr.matches(" *count\\( *\\) *" ) )
            {  functionName = "count";    // for backward compatibility
               fnargs = groups.get(0);    // just first group by column
               outTag = "ClusterCount";

            } else
            {  System.err.printf("Invalid aggregation function: %s\n", funcstr);
               exitWithHelp(options);
            }
         }

         try
         {  @SuppressWarnings("unchecked")
            Class<AggInterface> aFunctClass =
               (Class<AggInterface>) Class.forName(PACKAGEName + "." + functionName);
            Constructor<AggInterface> constr
                  = aFunctClass.getConstructor(String.class, String.class, String.class, String.class);

            ret[i] = constr.newInstance(outTag, resFormat, functionName, fnargs );

         } catch (InstantiationException e)
         {  System.err.println(e.getMessage());
            System.err.println("Function does not exist: " + functionName);
            exitWithHelp(options);

         } catch (IllegalAccessException e)
         {  System.err.println(e.getMessage());
            System.err.println("Function does not exist: " + functionName);
            exitWithHelp(options);

         } catch (ClassNotFoundException e)
         {  e.printStackTrace();
            System.err.println(e.getMessage());
            System.err.println("Function does not exist: " + functionName);
            exitWithHelp(options);

         } catch (NoSuchMethodException e)
         {  System.err.println("Function is not implemented correctly: " + functionName);
            System.err.println(e.getMessage());
            exitWithHelp(options);

         } catch (InvocationTargetException e)
         {  System.err.println(e.getTargetException().getMessage());
            exitWithHelp(options);

         } catch (Throwable e)
         {  System.err.println(e.getMessage());
            exitWithHelp(options);
         }
      }

      return ret;
   }

   /**
    *
    * Writes aggregated results to file.
    *
    */
   private static void writeGroup()
   {
      if (!currentGroupValue.isEmpty() && !oss.GetString().isEmpty())
      {
         oemolistream iss = new oemolistream();
         iss.SetFormat(OEFormat.SDF);

         oss.flush();
         // System.out.println(oss.GetString());
         if (!iss.openstring(oss.GetString()))
            oechem.OEThrow.Fatal("Unable to open input stream.");

         OEGraphMol mol2 = new OEGraphMol();
         int molCount = 0;
         while (oechem.OEReadMolecule(iss, mol2))
         {  if( outputMode == OUTPUTMODE.FIRST && molCount == 0 )
            {  outputMol(mol2, molCount);
               break;
            }
            else if( outputMode == OUTPUTMODE.ALL )
            {  outputMol(mol2, molCount);
            }
            else
            {  assert outputMode == OUTPUTMODE.LAST;
            }

            molCount++;
         }
         if( outputMode == OUTPUTMODE.LAST )
         {  outputMol(mol2, molCount-1);
         }

         mol2.delete();
         iss.close();
         iss.delete();
      }

      if (oss != null)
      {
         oss.close();
      }


   }

   private static void outputMol(OEGraphMol mol2, int molCount)
   {  for( AggInterface  fct : aAgg )
            oechem.OESetSDData(mol2, fct.getOutTagName(), fct.getResult(molCount));
         oechem.OEWriteMolecule(ofs, mol2);
   }

   /*
    * Process each molecule in a file by calling aggregator function.
    */
   private static void processRecord(OEGraphMol mol)
   {
      // Process current record
      String molGroupTagVal = getGroupTagValue(mol);
      // System.err.println(molGroupTagVal + " " + currentGroupValue);

      if (! molGroupTagVal.equals(currentGroupValue))
      {
         writeGroup();
         if (!oss.openstring())
            oechem.OEThrow.Fatal("Unable to open output stream.");
         for( AggInterface  fct : aAgg )
            fct.init();
         currentGroupValue = molGroupTagVal;
      }

      for( AggInterface  fct : aAgg )
         fct.process(mol);
      oechem.OEWriteMolecule(oss, mol);
   }

   private static String getGroupTagValue(OEGraphMol mol)
   {
      ArrayList<String> tagList = new ArrayList<String>();
      for (int i = 0; i < groups.size(); i++)
      {
         tagList.add(oechem.OEGetSDData(mol, groups.get(i)));
      }

      return tagList.toString();

   }

   public static void main(String[] args) throws IOException
   {  processOptions(args);

      ifs = new oemolistream();
      ofs = new oemolostream();
      oss = new oemolostream();
      oss.SetFormat(OEFormat.SDF);

      currentGroupValue = "";

      if (!ifs.open(inFile))
         oechem.OEThrow.Fatal("Unable to open " + inFile);
      if (!ofs.open(outFile))
         oechem.OEThrow.Fatal("Unable to create " + outFile);
      if (!oss.openstring())
         oechem.OEThrow.Fatal("Unable to open output stream.");

      OEGraphMol mol = new OEGraphMol();

      while (oechem.OEReadMolecule(ifs, mol))
      {  processRecord(mol);
      }
      writeGroup();

      ifs.close();
      ifs.delete();
      ofs.close();
      ofs.delete();
      oss.delete();
      mol.delete();
   }

   private static void exitWithHelp(Options options)
   {
      HelpFormatter formatter = new HelpFormatter();
      String head = "Note: input file must be ordered by the groupby tags";
      String foot = "\nFormat specs: siN (siginficant N digits) or rN (round to N digits)";
      formatter.printHelp("SDFAggregator", head, options, foot);
      //System.out.println("Format specs: siN (siginficant N digits) or rN (round to N digits)");
      System.exit(1);
   }
}
