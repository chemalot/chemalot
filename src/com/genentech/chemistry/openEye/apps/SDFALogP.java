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

import java.io.*;
import java.util.regex.Pattern;

import openeye.oechem.*;

import org.apache.commons.cli.*;

import com.aestel.utility.DataFormat;
import com.genentech.chemistry.openEye.aLogP.ALogPCalculator;
import com.genentech.oechem.tools.OETools;




/**

 * @author Alberto Gobbi/ 2012
 * Copyright 2012 Genentech
 */
public class SDFALogP
{  private static final String MY_NAME = "SDFALogP";
   private static final String ALOGP_TAG = "ALogP_GNE";
   
   private static final String OPT_INFILE    = "in";
   private static final String OPT_OUTFILE   = "out";
   private static final String OPT_SMARTS_FILE = "atomTypeFile";
   private static final String OPT_PRINT_COUNTS = "outputCounts";
   private static final String OPT_VALIDATE_ASSIGNMENT = "validateAssignment";
   private static final String OPT_SUPRESS_ZERO = "supressZeros";
   private static final String OPT_NEUTRALIZE = "neutralize";

   private final oemolothread outputOEThread;
   private final ALogPCalculator aLogPCalcualtor;
   private final boolean outputZero;
   private final boolean neutralize;
   
   private SDFALogP( String smartsFile, String outFile, boolean outputZero, 
            boolean neutralize, boolean validateAssignment )
   {  outputOEThread     = new oemolothread(outFile);
      this.outputZero = outputZero;
      this.neutralize = neutralize;
      aLogPCalcualtor = new ALogPCalculator( smartsFile, validateAssignment );
   }
   
   
   
   private void close()
   {  outputOEThread.close();
      aLogPCalcualtor.close();
   }
   
   

   private void run( String inFile, boolean outPutCounts )
   {  oemolithread ifs = new oemolithread(inFile);
      long start = System.currentTimeMillis();
      int iCounter = 0; //Structures in the SD file.

      OEMolBase mol = new OEGraphMol();
      OEMolBase compMol = new OEGraphMol();
      while( oechem.OEReadMolecule( ifs, mol ) )
      {  iCounter++;
    
         compMol.Clear();
         oechem.OEAddMols(compMol, mol);
         if( neutralize )
            OETools.neutralize(compMol);
            
         double alogp = aLogPCalcualtor.computeALogP(compMol);
         oechem.OEAddSDData(mol, ALOGP_TAG, DataFormat.formatNumber(alogp, "r2"));
         
         if( outPutCounts )
         {  int[] counts = aLogPCalcualtor.getAtomCounts();
            for(int i=1; i<counts.length; i++)
            {  if( outputZero || counts[i] > 0 )
                  oechem.OEAddSDData(mol, String.format("ALogP_GCount_%03d", i), 
                                          Integer.toString(counts[i])
                                          );
            }
         }
    
         oechem.OEWriteMolecule(outputOEThread, mol);
         
         //Output "." to show that the program is running.
         if( iCounter % 100 == 0 )
            System.err.print(".");
         if( iCounter % 4000 == 0 )
         {  System.err.printf( " %d %dsec\n",
                  iCounter, (System.currentTimeMillis()-start)/1000);
         }
      }

      mol.delete();
      compMol.delete();
      ifs.close();
      inFile = inFile.replaceAll( ".*" + Pattern.quote(File.separator), "" );
      System.err.printf( "%s: Read %d structures from %s. %d sec\n",
               MY_NAME, iCounter, inFile, (System.currentTimeMillis()-start)/1000 );
   }





   private static void exitWithHelp(Options options) 
   {  HelpFormatter formatter = new HelpFormatter();
      formatter.printHelp( MY_NAME, options );
      System.exit(1);
   }


   /**
    * @param args
    */
   public static void main( String...args ) throws IOException
   {  // create command line Options object
      Options options = new Options();
      Option opt = new Option( OPT_INFILE, true, 
               "input file oe-supported Use .sdf|.smi to specify the file type." );
      opt.setRequired( true );
      options.addOption( opt );

      opt = new Option( OPT_OUTFILE, true, 
               "output file oe-supported. Use .sdf|.smi to specify the file type." );
      opt.setRequired( true );
      options.addOption( opt );

      opt = new Option( OPT_SMARTS_FILE, true, 
               "Optional: to overwrite atom type definition file." );
      opt.setRequired( false );
      options.addOption( opt );

      opt = new Option( OPT_PRINT_COUNTS, false, 
               "If set the count of each atom type is added to the output file." );
      opt.setRequired( false );
      options.addOption( opt );

      opt = new Option( OPT_VALIDATE_ASSIGNMENT, false, 
               "Print warning if no atomtype matches an atom in a candidte molecule." );
      opt.setRequired( false );
      options.addOption( opt );

      opt = new Option( OPT_SUPRESS_ZERO, false, 
               "If given atom type counts with count=0 will not be added." );
      opt.setRequired( false );
      options.addOption( opt );

      opt = new Option( OPT_NEUTRALIZE, true, 
               "y|n to neutralize molecule if possible (default=y)" );
      opt.setRequired( false );
      options.addOption( opt );

      CommandLineParser parser = new PosixParser();
      CommandLine cmd = null;
      try
      {  cmd = parser.parse( options, args );
      } catch( Exception e )
      {  System.err.println( e.getMessage() );
         exitWithHelp( options );
      }
      args = cmd.getArgs();

      if( cmd.hasOption( "d" ) )
      {  System.err.println( "Start debugger and press return:" );
         new BufferedReader( new InputStreamReader( System.in ) ).readLine();
      }

      String inFile = cmd.getOptionValue( OPT_INFILE );
      String outFile = cmd.getOptionValue( OPT_OUTFILE );
      String smartsFile = cmd.getOptionValue( OPT_SMARTS_FILE );
      boolean outputCount = cmd.hasOption(OPT_PRINT_COUNTS);
      boolean outputZero = ! cmd.hasOption(OPT_SUPRESS_ZERO);
      boolean neutralize = ! "n".equalsIgnoreCase(cmd.getOptionValue(OPT_NEUTRALIZE));
      boolean ValidateAssignment = cmd.hasOption(OPT_VALIDATE_ASSIGNMENT);
      
      SDFALogP sdfALogP= new SDFALogP( smartsFile, outFile, outputZero, neutralize, ValidateAssignment  );
      
      sdfALogP.run( inFile, outputCount );
      sdfALogP.close();
   }
}
