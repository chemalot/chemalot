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
package com.genentech.chemistry.tool;


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
import com.aestel.utility.Message;
import com.aestel.utility.exception.IncorrectInputException;



public class SDFSelectivityCalculator
{  private static final String MY_NAME = "SDFSelectivityCalculator";
   private static final String OPT_INFILE         = "in";
   private static final String OPT_OUTFILE        = "out";
   private static final String OPT_NUMERATOR      = "numerator";
   private static final String OPT_NUMERATOR_OP   = "numeratorOp";
   private static final String OPT_DENOMINATOR    = "denominator";
   private static final String OPT_DENOMINATOR_OP = "denominatorOp";
   private static final String OPT_SELECTIVITY    = "selectivity";
   private static final String OPT_SELECTIVITY_OP = "selectivityOp";
   private static final String OPT_OUTPUT_MODE    = "outputMode";

   public static enum OutputMode
   {  SEPARATE( "separate" ),
      COMBINE( "combine" ),
      COMBINE_OP( "combine+op" );

      private final String value;
      OutputMode( String value )
      {  this.value = value;
      }

      @Override
      public String toString()
      {  return value; }

      public static OutputMode toOutputMode( String str )
      {  OutputMode[] modes = OutputMode.values();
         for( OutputMode mode : modes )
         {  if( mode.toString().equalsIgnoreCase( str ) )
               return mode;
         }
         return null;
      }
   }
   
   
   private final oemolothread outputOEThread;
   private String numeratorTag;
   private String numeratorOPTag;
   private String denominatorTag;
   private String denominatorOPTag;
   private String selectivityTag;
   private String selectivityOPTag;
   private String messageTag;

   private OutputMode outputMode;
   private String numeratorInput;
   private String numeratorOPInput;
   private String denominatorInput;
   private String denominatorOPInput;
   private String selectivity;
   private String selectivityOP;
   private String messageText;

   
   private SDFSelectivityCalculator( String outFile ) 
   {  this.outputOEThread = new oemolothread(outFile);  }
   
   
   private void setParameters( String outputMode,
            String numeratorTag, String numeratorOPTag, 
            String denominatorTag, String denominatorOPTag, 
            String selectivityTag, String selectivityOPTag)
   throws IncorrectInputException
   {  this.outputMode = OutputMode.COMBINE_OP;
      if( outputMode != null && outputMode.length() > 0 )
      {  this.outputMode = OutputMode.toOutputMode( outputMode );
         if( this.outputMode == null )
         {  String s = outputMode + " is not valid";
            Message m = new Message(s, Message.Level.ERROR, null );
            throw new IncorrectInputException( m );
         }
      }
      this.numeratorTag   = numeratorTag;
      this.denominatorTag = denominatorTag;
      
      if( numeratorOPTag == null || numeratorOPTag.length() == 0 )
         this.numeratorOPTag = null;
      else
         this.numeratorOPTag = numeratorOPTag;

      if( denominatorOPTag == null || denominatorOPTag.length() == 0 )
         this.denominatorOPTag = null;
      else
         this.denominatorOPTag = denominatorOPTag;
   
      if( selectivityTag == null || selectivityTag.length() == 0 )
         this.selectivityTag = denominatorTag + " Selectivity";
      else
         this.selectivityTag = selectivityTag;
      
      if( selectivityOPTag == null || selectivityOPTag.length() == 0 )
         this.selectivityOPTag = this.selectivityTag + " OP";
      else
         this.selectivityOPTag = selectivityOPTag;
      this.messageTag = this.selectivityTag + " Error";
   }
   
   
   private void close()
   {  outputOEThread.close();  }


   private void run( String inFile )
   {  oemolithread ifs = new oemolithread( inFile );
      long start = System.currentTimeMillis();
      int iCounter = 0; //Structures in the SD file.

      OEMolBase mol = new OEGraphMol();
      while( oechem.OEReadMolecule( ifs, mol ) )
      {  iCounter++;
         if( readInput( mol ) )
         {  calculate();
         
            if( outputMode == OutputMode.SEPARATE )
            {  oechem.OESetSDData( mol, selectivityOPTag, selectivityOP );
               oechem.OESetSDData( mol, selectivityTag, selectivity );
            } else 
            {  oechem.OESetSDData( mol, selectivityTag, selectivityOP + selectivity );
               if( outputMode == OutputMode.COMBINE_OP )
                  oechem.OESetSDData( mol, selectivityOPTag, selectivityOP );
            }
         } else
         {  if( messageText != null && messageText.length() > 0 )
               oechem.OESetSDData( mol, messageTag, messageText );
         }
         oechem.OEWriteMolecule( outputOEThread, mol );
      }

      //Output "." to show that the program is running.
      if( iCounter % 100 == 0 )
         System.err.print(".");
      if( iCounter % 4000 == 0 )
      {  System.err.printf( " %d %dsec\n",
               iCounter, ( System.currentTimeMillis()-start )/1000 );
      }
      mol.delete();
      ifs.close();
      inFile = inFile.replaceAll( ".*" + Pattern.quote( File.separator ), "" );
      System.err.printf( "%s: Read %d structures from %s. %d sec\n",
               MY_NAME, iCounter, inFile, ( System.currentTimeMillis()-start )/1000 );
   }
   
   
   private boolean readInput( OEMolBase mol )
   {
      numeratorInput     = oechem.OEGetSDData( mol, numeratorTag );
      denominatorInput   = oechem.OEGetSDData( mol, denominatorTag );
      if( numeratorInput == null || numeratorInput.length() == 0 || 
          denominatorInput == null || denominatorInput.length() == 0 )
         return false; 
      
      numeratorOPInput = "";
      if( numeratorOPTag != null )
         numeratorOPInput   = oechem.OEGetSDData( mol, numeratorOPTag );
      
      denominatorOPInput = "";
      if( denominatorOPTag != null )
      denominatorOPInput = oechem.OEGetSDData( mol, denominatorOPTag );
      
      if( numeratorOPInput.length() == 0 )
      {  if( numeratorInput.startsWith( "<" ) || numeratorInput.startsWith( ">" ) )
         {  numeratorOPInput = numeratorInput.substring( 0, 1 );
            numeratorInput   = numeratorInput.substring( 1 );
         }
      }
      if( denominatorOPInput.length() == 0 )
      {  if( denominatorInput.startsWith( "<" ) || denominatorInput.startsWith( ">" ) )
         {  denominatorOPInput = denominatorInput.substring( 0, 1 );
            denominatorInput = denominatorInput.substring( 1 );
         }
      }
      
      StringBuffer errText = new StringBuffer();
      if( !DataFormat.isNumeric( numeratorInput ) || 
          !DataFormat.isNumeric( denominatorInput ) )
      {   errText.append( "Numerator or denominator is not a number; ");
      }
      if( numeratorOPInput.length() > 0 && !numeratorOPInput.matches( "[<>]" ) )
      {   errText.append( String.format("Invalid numerator operators (%s); ", numeratorOPInput));
      }
      if( denominatorOPInput.length() > 0 && !denominatorOPInput.matches( "[<>]" ) )
      {   errText.append( String.format("Invalid denominator operators (%s); ", denominatorOPInput));
      }
      if( errText.length() > 0 )
      {  messageText = errText.substring(0, errText.length()-2);
         return false;
      }
      return true;
   }

   
   private void calculate()
   {  float  numerator;
      float  denominator;
      
      if( numeratorInput == null || numeratorInput.length() == 0  || 
          denominatorInput == null || denominatorInput.length() == 0 )
      {  selectivityOP = "";
         selectivity = "";
         return;
      }
      
      selectivityOP = numeratorOPInput;
      if( numeratorOPInput.length() == 0 && denominatorOPInput.length() == 0 )
      {  selectivityOP = "";
      } else if( numeratorOPInput.length() == 0 && denominatorOPInput.length() > 0 )
      {  if( denominatorOPInput.equals( ">" ) )
         selectivityOP = "<";
         else
            selectivityOP = ">";
      } else if( numeratorOPInput.length() > 0 && denominatorOPInput.length() == 0 )
      {  selectivityOP = numeratorOPInput;
      } else if( numeratorOPInput.equals( denominatorOPInput ) )
      {  selectivityOP =  "u";
      }
      
      numerator = Float.parseFloat( numeratorInput );
      denominator = Float.parseFloat( denominatorInput );
      selectivity = DataFormat.formatNumber( numerator/denominator, "si2" );
   }
   
   
   private static void exitWithHelp( Options options ) 
   {  String h = "Calculate selectivity with consideration of the operators.";
      String s = MY_NAME + " -in fn -out fn  -outputMode combine+op" +
                 " -numerator fieldName -numeratorOp fieldName" +
                 " -denominator fieldName -denominatorOp fieldName" +
                 " -selectivity fieldName -selectivityOp fieldName";
      HelpFormatter formatter = new HelpFormatter();
      formatter.printHelp( s, h, options, null );
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
      
      opt = new Option( OPT_NUMERATOR, true, 
               "Name of the field containing the numerator value." );
      opt.setRequired( true );
      options.addOption( opt );
      
      opt = new Option( OPT_NUMERATOR_OP, true, 
               "Name of the field containing the numerator operator." +
               " If this option is not set it is assumed that the operator" +
               " is included in the " + OPT_NUMERATOR + " field." );
      opt.setRequired( false );
      options.addOption( opt );
      
      opt = new Option( OPT_DENOMINATOR, true, 
               "Name of the field containing the denominator value." );
      opt.setRequired( true );
      options.addOption( opt );
      
      opt = new Option( OPT_DENOMINATOR_OP, true, 
               "Name of the field containing the denominator operator." +
               " If this option is not set it is assumed that the operator" +
               " is included in the " + OPT_DENOMINATOR + " field." );
      opt.setRequired( false );
      options.addOption( opt );
      
      opt = new Option( OPT_SELECTIVITY, true, 
               "Name of the field containing the output selectivity value." +
               " If this option is not set, the field name of the denominator" +
               " is used." );
      opt.setRequired( false );
      options.addOption( opt );
      
      opt = new Option( OPT_SELECTIVITY_OP, true, 
               "Name of the field containing the output selectivity operator." +
               " If this option is not set it is assumed that the operator" +
               " is included in the " + OPT_SELECTIVITY + " field." );
      opt.setRequired( false );
      options.addOption( opt );
      
      opt = new Option( OPT_OUTPUT_MODE, true, 
               "Valid: separate, combine, combine+op;" +
               " If not specified, default is combine+op." +
               " (separate=operator and value in separate fields;" +
               " combine=operator and value in one field;" + 
               " combine+op=one field with the operator and value and" +
               " plus one field with only the operator)");
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
   
      String inFile           = cmd.getOptionValue( OPT_INFILE );
      String outFile          = cmd.getOptionValue( OPT_OUTFILE );
      String outputMode       = cmd.getOptionValue( OPT_OUTPUT_MODE );
      String numeratorTag     = cmd.getOptionValue( OPT_NUMERATOR );
      String numeratorOPTag   = cmd.getOptionValue( OPT_NUMERATOR_OP );
      String denominatorTag   = cmd.getOptionValue( OPT_DENOMINATOR );
      String denominatorOPTag = cmd.getOptionValue( OPT_DENOMINATOR_OP );
      String selectivityTag   = cmd.getOptionValue( OPT_SELECTIVITY );
      String selectivityOPTag = cmd.getOptionValue( OPT_SELECTIVITY_OP );
      
      SDFSelectivityCalculator calculator = new SDFSelectivityCalculator( outFile );
      try
      {  calculator.setParameters( 
               outputMode, numeratorTag, numeratorOPTag, denominatorTag, 
               denominatorOPTag, selectivityTag, selectivityOPTag );
         calculator.run( inFile );
      } catch( IncorrectInputException e )
      {  System.err.println( e.getMessage() );
         exitWithHelp( options );
      } finally
      {  calculator.close();
      }
   }
}
