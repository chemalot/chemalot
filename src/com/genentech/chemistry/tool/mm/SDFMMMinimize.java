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
package com.genentech.chemistry.tool.mm;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.PosixParser;


/**

 * @author Ben Sellers / Alberto Gobbi / 2014
 * Copyright 2014 Genentech
 * 
 * Class to minimize input molecules with a selected molecular
 *   mechanics force field.
 *  
 */
public class SDFMMMinimize
{
   private static final String MY_NAME = "sdfMMMinimize";
   
   private static final String OPT_INFILE         = "in";
   private static final String OPT_OUTFILE        = "out";
   private static final String OPT_PROGRAM        = "program";
   private static final String OPT_FORCEFIELD     = "forcefield";
   private static final String OPT_SOLVENT        = "solvent";
   private static final String OPT_FIXED_ATOM_TAG = "fixedAtomTag";
   private static final String OPT_FIX_TORSION    = "fixTorsion";   
   private static final String OPT_WORKING_DIR    = "workDir";
   
   private MMMinMethod minMethod = null;
      
   private final static Map<String,MMMinMethod> minMethodNameObjectMap;
   static
   {
      //
      // We are creating all lightweight minimizers statically so we have access to 
      //  the available command line options for each which are defined at 
      //  object creation time.
      //
      minMethodNameObjectMap = new HashMap<String, MMMinMethod>();
      minMethodNameObjectMap.put("SZYBKI",    new SZYBKIMinMethod());
      minMethodNameObjectMap.put("MOE",       new MOEMinMethod());
      minMethodNameObjectMap.put("MACROMODEL",new MacromodelMinMethod());

      /*************************************************
       * ADD NEW MM METHOD HERE: list of known programs
       *************************************************/
   }
   
   /*
    * Create a minimization method instance matching the requested program and forcefield names
    * @param requestedProgram the name of the program to use
    * @param requestedFF the name of the forcefield to use
    */
   private static MMMinMethod initializeMinimizeMethod (String requestedProgram, String requestedFF, String requestedSolvent)
   {
      MMMinMethod minMethod = null;
      
      for (String knownProgram : minMethodNameObjectMap.keySet())
      {
         if (knownProgram.equalsIgnoreCase(requestedProgram))
         {                        
            minMethod = getMinMethodFromName(requestedProgram);
            
            if (requestedFF != null)
            {    minMethod.setForceField(requestedFF);    }
            
            if (requestedSolvent != null)
            {    minMethod.setSolvent(requestedSolvent);    }
            break;
         }
      }
      if (minMethod == null)
      { throw new Error ("Minimization method not available." + requestedProgram);    }
      return minMethod;
   }
   
   /*
    * Get the minimization method object from a program name
    * @param requestedProgram the requested program name
    * @return a minimization instance wrapping that program
    */
   private static MMMinMethod getMinMethodFromName(String requestedProgram)
   {
      if (minMethodNameObjectMap.containsKey(requestedProgram))
      {
         return minMethodNameObjectMap.get(requestedProgram);
      }
      { 
         throw new Error ("Minimization method not available." + requestedProgram);    
      }
   }

   /*
    * Default constructior
    */
   public SDFMMMinimize()
   {     
   }
   
   /*
    * Set which program, forcefield and solvent are to be used
    * @param program name
    * @param forcefield option name
    * @param solvent option name
    */
   public void setMethod (String program, String forcefield, String solvent)
   {
      this.minMethod = initializeMinimizeMethod(program, forcefield, solvent);
   }

   /**
    * Main function for running on the command line
    * @param args
    */
   public static void main( String...args ) throws IOException
   {  
      // Get the available options from the programs
      Map<String,List<String>> allowedProgramsAndForceFields = getAllowedProgramsAndForceFields(); 
      Map<String,List<String>> allowedProgramsAndSolvents    = getAllowedProgramsAndSolvents();      
      
      // create command line Options object
      Options options = new Options();
      Option opt = new Option( OPT_INFILE, true, 
               "input file oe-supported Use .sdf|.smi to specify the file type." );
      opt.setRequired( true );
      options.addOption( opt );

      opt = new Option( OPT_OUTFILE, true, 
               "output file oe-supported. Use .sdf|.smi to specify the file type." );
      opt.setRequired( true );
      options.addOption( opt );
      
      StringBuilder programOptions = new StringBuilder("Program to use for minimization.  Choices are\n");
      
      for (String program : allowedProgramsAndForceFields.keySet())
      {
         programOptions.append("\n***   -program " + program + "   ***\n");
         String forcefields = "";
         for (String option : allowedProgramsAndForceFields.get(program))
         {
            forcefields += option + " " ;           
         }
         programOptions.append("-forcefield " + forcefields + "\n");
         
         String solvents = "";
         for (String option : allowedProgramsAndSolvents.get(program))
         {
            solvents += option + " ";
         }
         programOptions.append("-solvent " + solvents + "\n");
      }
            
      opt = new Option( OPT_PROGRAM, true, 
               programOptions.toString() );
      opt.setRequired( true );
      options.addOption( opt );
      
      opt = new Option( OPT_FORCEFIELD, true, 
               "Forcefield options.  See -program for choices" );
      opt.setRequired( false );
      options.addOption( opt );
      
      opt = new Option( OPT_SOLVENT, true, 
               "Solvent options.  See -program for choices" );
      opt.setRequired( false );
      options.addOption( opt );
      
      opt = new Option( OPT_FIXED_ATOM_TAG, true, 
               "SD tag name which contains the atom numbers to be held fixed." );
      opt.setRequired( false );
      options.addOption( opt );    
      
      opt = new Option( OPT_FIX_TORSION, true, 
               "true/false. if true, the atoms in fixedAtomTag contains 4 indices of atoms defining a torsion angle to be held fixed" );
      opt.setRequired( false );
      options.addOption( opt );    
      
      opt = new Option( OPT_WORKING_DIR, true, 
               "Working directory to put files." );
      opt.setRequired( false );
      options.addOption( opt );    
      
      CommandLineParser parser = new PosixParser();
      CommandLine cmd = null;
      try
      {  
         cmd = parser.parse( options, args );
      } 
      catch( Exception e )
      {  
         System.err.println( e.getMessage() );
         exitWithHelp( options );
      }
      args = cmd.getArgs();

      if (args.length != 0)
      {  System.err.println("Unknown arguments" + args);
         exitWithHelp(options);
      }

      if (cmd.hasOption("d"))
      {  System.err.println("Start debugger and press return:");
         new BufferedReader(new InputStreamReader(System.in)).readLine();
      }

      String inFile       = cmd.getOptionValue( OPT_INFILE );
      String outFile      = cmd.getOptionValue( OPT_OUTFILE );
      String fixedAtomTag = cmd.getOptionValue( OPT_FIXED_ATOM_TAG );
      boolean fixTorsion  = ( cmd.getOptionValue( OPT_FIX_TORSION ) != null 
                           && cmd.getOptionValue(OPT_FIX_TORSION).equalsIgnoreCase("true"));
      String programName  = cmd.getOptionValue( OPT_PROGRAM ); 
      String forcefield   = cmd.getOptionValue( OPT_FORCEFIELD );  
      String solvent      = cmd.getOptionValue( OPT_SOLVENT );    
      String workDir      = cmd.getOptionValue( OPT_WORKING_DIR );           
      
      if (workDir == null || workDir.trim().length() == 0) 
         workDir = ".";
      
      // Create a minimizer 
      SDFMMMinimize minimizer = new SDFMMMinimize ();
      minimizer.setMethod(programName, forcefield, solvent);
      minimizer.run (inFile, outFile, fixedAtomTag, fixTorsion, workDir);
      minimizer.close();
      System.err.println("Minimization complete.");
   }

   /*
    * Internal method to get all the allowed forcefield names
    * @return a map of program->list of forcefield names
    */
   private static Map<String, List<String>> getAllowedProgramsAndForceFields()
   {
      Map<String, List<String>> programOptionMap = new HashMap<String, List<String>>();
      
      for (String program :  minMethodNameObjectMap.keySet())
      {
         MMMinMethod minMethod = getMinMethodFromName(program);
         if (minMethod != null)
         {
            List<String> forcefields =  Arrays.asList(minMethod.getAvailableForceFieldNames());
            programOptionMap.put(minMethod.getMethodName(), forcefields);
         }
      }
      return programOptionMap;
   }
   
   /*
    * Internal method to get all the allowed solvent names
    * @return a map of program->list of solvent names
    */
   private static Map<String, List<String>> getAllowedProgramsAndSolvents()
   {
      Map<String, List<String>> programOptionMap = new HashMap<String, List<String>>();
      
      for (String program :  minMethodNameObjectMap.keySet())
      {
         MMMinMethod minMethod = getMinMethodFromName(program);
         if (minMethod != null)
         {
            List<String> solvents =  Arrays.asList(minMethod.getAvailableSolventNames());         
            programOptionMap.put(minMethod.getMethodName(), solvents);
         }
      }
      return programOptionMap;
   }

   /*
    * Print the help message
    * @param options
    */
   private static void exitWithHelp(Options options)
   {
      HelpFormatter formatter = new HelpFormatter();
      String head = "Minimizes input molecules with a specified program and force-field.  Select atoms can be fixed.";
      formatter.printHelp( MY_NAME, head, options, "", true );
      System.exit(1);
   }

   /*
    * Run the minimization
    * @param inFile the input filename
    * @param outFile the output filename
    * @param fixedAtomTag the sdf tag name that contains the atom indices of torsion to be held fixed
    * @param fixTorsion fixedAtomTag lists atom indices which define a four-atom torsion angle to be held fixed if true
    * @param workDir the working directory full path
    */
   private void run(String inFile, String outFile, String fixedAtomTag, boolean fixTorsion, String workDir)
   {       
      minMethod.execute(inFile, outFile, fixedAtomTag, fixTorsion, workDir);
   }

   
   /*
    * Close the session
    */
   public void close() 
   { 
      //outputOEThread.close();
   }
}
