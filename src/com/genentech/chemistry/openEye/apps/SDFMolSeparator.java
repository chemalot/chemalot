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

import com.genentech.chemistry.openEye.MolSeparator;

/**
 * Separate components in a molfile into individual entities which could be
 * as sdf record, smiles, or whichever file format supported by OEChem Toolkits.
 *
 * @author Man-Ling Lee / December 20, 2012
 * Copyright 2012-2013 Genentech
 */
public class SDFMolSeparator
{
   private static final String MY_NAME = "sdfMolSeparator";
   private static final String OPT_INFILE    = "in";
   private static final String OPT_OUTFILE   = "out";


   private final oemolothread outputOEThread;
   private final MolSeparator molSeparator;


   private SDFMolSeparator( String outputFile )
   {
      outputOEThread = new oemolothread( outputFile );
      molSeparator = new MolSeparator();
   }


   private void close()
   {  outputOEThread.close();
      outputOEThread.delete();
      molSeparator.clear();
   }


   private void run( String inFile )
   {
      oemolithread ifs = new oemolithread( inFile );
      long start = System.currentTimeMillis();
      int iCounter = 0; //Structures in the input SD file
      int oCounter = 0; //Structure in the output SD file

      OEMolBase inMol = new OEGraphMol();
      while( oechem.OEReadMolecule( ifs, inMol ) )
      {
         iCounter++;
         int nComponents = molSeparator.separate( inMol );
         for( int i=0; i<nComponents; i++ )
         {
            OEMolBase mol = molSeparator.getMol( i );
            oechem.OEWriteMolecule( outputOEThread, mol );
            oCounter++;
         }
         molSeparator.clear();

         //Output "." to show that the program is running.
         if( iCounter % 100 == 0 )
            System.err.print(".");
         if( iCounter % 6000 == 0 )
         {  System.err.printf( " %d %dsec\n",
                  iCounter, (System.currentTimeMillis()-start)/1000);
         }
      }

      inMol.delete();
      ifs.close();
      ifs.delete();
      inFile = inFile.replaceAll( ".*" + Pattern.quote(File.separator), "" );
      System.err.printf(
               "%s: Read %d molecules from %s and output %d molecules. %d sec\n",
               MY_NAME, iCounter, inFile, oCounter,
               (System.currentTimeMillis()-start)/1000 );
   }



   private static void exitWithHelp( Options options )
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
      SDFMolSeparator separator = new SDFMolSeparator( outFile );

      try
      {  separator.run( inFile );

      } catch( IndexOutOfBoundsException iie )
      {
         System.err.println( iie.toString() );
         exitWithHelp( options );
      } finally
      {  separator.close();
      }

   }
}
