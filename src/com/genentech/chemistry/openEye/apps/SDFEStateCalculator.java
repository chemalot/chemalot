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

import com.aestel.utility.DataFormat;
import com.aestel.utility.NameValuePair;
import com.genentech.chemistry.openEye.EState.EStateCalculator;



/**
 * Parsing SD file or file with format accepted by OEChem tool kit.
 *
 * @author Man-Ling Lee / August 04, 2012
 * Copyright 2012-2015 Genentech
 */
public class SDFEStateCalculator
{  private static final String MY_NAME =  SDFEStateCalculator.class.getSimpleName();
   private static final String OPT_INFILE           = "in";
   private static final String OPT_OUTFILE          = "out";
   private static final String OPT_ESTATE_COUNT     = "es_count";
   private static final String OPT_ESTATE_SUM       = "es_sum";
   private static final String OPT_ESTATE_SYMBOL    = "es_symbol";
   private static final String OPT_UNASSIGNED_ATOMS = "unassigned_atoms";
   private static final String OPT_UNASSIGNED_COUNT = "unassigned_count";
   private static final String OPT_ESTATE_INDICE    = "es_indice";
   private static final String OPT_SMARTS           = "smarts";
   private static final String OPT_PRINT_DETAILS    = "print_details";
   
   private static final String TAG_ES_COUNT         = "ES_Count";
   private static final String TAG_ES_SUM           = "ES_Sum";
   private static final String TAG_ES_SYMBOL        = "ES_Symbol";
   private static final String TAG_ESTATE_PREFIX    = "EState";
   private static final String TAG_UNASSIGNED_COUNT = "ES_Unassigned_Count";
   private static final String TAG_UNASSIGNED_ATOMS = "ES_Unassigned_Atoms";
   
   private static void exitWithHelp( Options options ) 
   {  HelpFormatter formatter = new HelpFormatter();
      formatter.printHelp( MY_NAME, options );
      System.exit(1);
   }

   

   private final oemolothread outputOEThread;
   private final EStateCalculator esCalculator;
   private boolean outputESCount;
   private boolean outputESSum;
   private boolean outputESSymbol;
   private boolean outputUnassignedCount;
   private boolean outputUnassignedAtoms;
   private boolean outputESIndex;
   private String  smarts;
   private boolean printDetails;
   
   
   private SDFEStateCalculator( String outFile )
   {  outputOEThread = new oemolothread( outFile );
      esCalculator   = new EStateCalculator();
   }
   
   
   private void prepare( boolean outputESCount, boolean outputESSum, 
                         boolean outputESSymbol, boolean outputUnknownCount, 
                         boolean outputUnknownAtoms, boolean outputESIndex, 
                         String smarts, boolean printDetails )
   {  this.outputESCount  = outputESCount;
      this.outputESSum    = outputESSum;
      this.outputESSymbol = outputESSymbol;
      this.outputUnassignedCount = outputUnknownCount;
      this.outputUnassignedAtoms = outputUnknownAtoms;
      this.outputESIndex  = outputESIndex;
      this.smarts         = smarts;
      this.printDetails   = printDetails;
   }
   
   
   private void run( String inFile )
   {  oemolithread ifs = new oemolithread( inFile );
      long start = System.currentTimeMillis();
      int iCounter = 0; //Structures in the SD file.

      OEMolBase mol = new OEGraphMol();
      while( oechem.OEReadMolecule( ifs, mol ) )
      {  iCounter++;
         esCalculator.compute( mol, outputESSum, printDetails );
         
         if( outputESIndex )
            oechem.OESetSDData( mol, TAG_ESTATE_PREFIX + "_INDICE", 
                                esCalculator.getEStateIndexSummary() );
         if( smarts != null)
         {  NameValuePair<String,String>[] nvPairs = esCalculator.getEStateOf( mol, smarts );
            addOutputs( mol, TAG_ESTATE_PREFIX, nvPairs );
         }
         if( outputESCount )
            addOutputs( mol, TAG_ES_COUNT, 
                        esCalculator.getEStateCounts() );
         if( outputESSum )
           addOutputs( mol, TAG_ES_SUM, 
                       esCalculator.getEStateSums() );
         if( outputESSymbol )
            addOutputs( mol, TAG_ES_SYMBOL, 
                        esCalculator.getEStateAtomGroupSymbols() );
         if( outputUnassignedCount )
            oechem.OESetSDData( mol, TAG_UNASSIGNED_COUNT, 
                                String.valueOf( esCalculator.getUnknownCount() ) );
         if( outputUnassignedAtoms )
            oechem.OESetSDData( mol, TAG_UNASSIGNED_ATOMS, 
                                esCalculator.getUnknownAtoms() );
         oechem.OEWriteMolecule( outputOEThread, mol );
      }

      //Output "." to show that the program is running.
      if( iCounter % 100 == 0 )
         System.err.print(".");
      if( iCounter % 4000 == 0 )
      {  System.err.printf( " %d %dsec\n",
               iCounter, (System.currentTimeMillis()-start)/1000);
      }
      mol.delete();
      ifs.close();
      ifs.delete();
      inFile = inFile.replaceAll( ".*" + Pattern.quote(File.separator), "" );
      System.err.printf( "%s: Read %d structures from %s. %d sec\n",
               MY_NAME, iCounter, inFile, (System.currentTimeMillis()-start)/1000 );
   }
   
   
   private void addOutputs( OEMolBase mol, String tagPrefix, int[] values )
   {  for( int i=0; i<values.length; i++ ) 
      {  String tag = tagPrefix + "_" + formatNumber( i+1 );
         oechem.OESetSDData( mol, tag, String.valueOf( values[i] ) );
      }
   }
   
   private void addOutputs( OEMolBase mol, String tagPrefix, float[] values )
   {  for( int i=0; i<values.length; i++ ) 
      {  String tag = tagPrefix + "_" + formatNumber( i+1 );
         String val = DataFormat.formatNumber(values[i], "r3" );
         oechem.OESetSDData( mol, tag, val );
      }
   }
   
   private void addOutputs( OEMolBase mol, String tagPrefix, String[] values )
   {  for( int i=0; i<values.length; i++ ) 
      {  String tag = tagPrefix + "_" + formatNumber( i+1 );
         oechem.OESetSDData( mol, tag, values[i] );
      }
   }
   
   static private void addOutputs( OEMolBase mol, String tagPrefix, NameValuePair<String,String>[] nvPairs )
   {  for( int i=0; i<nvPairs.length; i++ ) 
      {  String tag = tagPrefix + "_" + nvPairs[i].getName();
         oechem.OESetSDData( mol, tag, nvPairs[i].getValue() );
      }
   }
   
   @SuppressWarnings("static-method")
   private String formatNumber ( int number )
   {  if( number < 10 )
         return "0" + String.valueOf( number );
      return String.valueOf( number );
   }
   
   
   private void close()
   {  outputOEThread.close();
      outputOEThread.delete();
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

      opt = new Option( OPT_ESTATE_COUNT, false, 
               "Output the the counts (occurrence) for each E-state atom group." +
               " E-state counts will be output if no output option is specified" );
      opt.setRequired( false );
      options.addOption( opt );

      opt = new Option( OPT_ESTATE_SUM, false, 
               "Output the sum of the E-state indices for each E-state atom group" +
               " in addition to the output of the" + 
               " E-state atom groups." );
      opt.setRequired( false );
      options.addOption( opt );

      opt = new Option( OPT_ESTATE_SYMBOL, false, 
               "Output the E-state atom group symbol" );
      opt.setRequired( false );
      options.addOption( opt );

      opt = new Option( OPT_UNASSIGNED_ATOMS, false, 
               "Output atoms in the given molecules that" +
               " could not be assigned to an E-state atom group" );
      opt.setRequired( false );
      options.addOption( opt );

      opt = new Option( OPT_UNASSIGNED_COUNT, false, 
               "Output the number of atoms in the given molecules that" +
               " could not be assigned to an E-state atom group" );
      opt.setRequired( false );
      options.addOption( opt );

      opt = new Option( OPT_ESTATE_INDICE, false, 
               "Output the E-state indice of all atoms in the given molecule" + 
               " in one field." );
      opt.setRequired( false );
      options.addOption( opt );

      opt = new Option( OPT_SMARTS, true, 
               "Description of the group of interest. If specified, the" + 
               " E-state indice of atoms matching the SMARTS are output" +
               " in separated fields." );
      opt.setRequired( false );
      options.addOption( opt );

      opt = new Option( OPT_PRINT_DETAILS, false, 
               "Output the details of the calculation to stderr" );
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

      String inFile  = cmd.getOptionValue( OPT_INFILE );
      String outFile = cmd.getOptionValue( OPT_OUTFILE );
      String smarts  = cmd.getOptionValue( OPT_SMARTS );
      boolean outputESIndex  = cmd.hasOption( OPT_ESTATE_INDICE );
      boolean outputESCount  = cmd.hasOption( OPT_ESTATE_COUNT );
      boolean outputESSum    = cmd.hasOption( OPT_ESTATE_SUM );
      boolean outputESSymbol = cmd.hasOption( OPT_ESTATE_SYMBOL );
      boolean outputUnkCount = cmd.hasOption( OPT_UNASSIGNED_COUNT );
      boolean outputUnkAtoms = cmd.hasOption( OPT_UNASSIGNED_ATOMS );
      boolean printDetails   = cmd.hasOption( OPT_PRINT_DETAILS );
      SDFEStateCalculator calculator = new SDFEStateCalculator( outFile );
      
      if( !outputESCount && !outputESSum && !outputUnkCount && !outputESSymbol
       && !outputESIndex && ( smarts == null || smarts.length() == 0 ) )
         outputESCount = true;
      
      try
      {  calculator.prepare( outputESCount, outputESSum, outputESSymbol, 
                             outputUnkCount, outputUnkAtoms, outputESIndex, 
                             smarts, printDetails );
         calculator.run( inFile );
      } finally
      {  calculator.close();
      }
   }
   
}
