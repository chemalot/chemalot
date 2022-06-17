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

import com.genentech.chemistry.openEye.RingSystemExtractor;


/**
 * Identify the rings systems of the input structure and append fields with
 * largest ring system, largest ring, and the basic ring systems to the SD file.
 *
 * @author Man-Ling Lee / August 2, 2010
 * Copyright 2010 Genentech
 */
public class SDFRingSystemExtractor
{  private static final String MY_NAME = "SDFRingSystemExtraction";
   private static final String OPT_INFILE    = "in";
   private static final String OPT_OUTFILE   = "out";
   private static final String OPT_REPLACE_OUTPUT_WITH_LARGEST_RING_SYSTEM = "replaceOutputWithLargestRingSystem";
   private static final String OPT_REPLACE_OUTPUT_WITH_LARGEST_RING        = "replaceOutputWithLargestRing";
   private static final String OPT_KEEP_ONE_SIDE_CHAIN_ATOM                = "keepOneSideChainAtom";
   
   private static boolean replaceOutputWithLargestRingSystem = false;
   private static boolean replaceOutputWithLargestRing       = false;
   private static boolean keepOneSideChainAtom               = false;
   
   private final oemolothread outputOEThread;


   private SDFRingSystemExtractor( String outFile )
   {  outputOEThread = new oemolothread(outFile);
   }


   private void close()
   {  outputOEThread.close();
   }


   private void run( String inFile )
   {  oemolithread ifs = new oemolithread(inFile);
      long start = System.currentTimeMillis();
      int iCounter = 0; //Structures in the SD file.

      OEMolBase mol = new OEGraphMol();
      RingSystemExtractor rsExtractor = new RingSystemExtractor(keepOneSideChainAtom);
      while( oechem.OEReadMolecule( ifs, mol ) )
      {  iCounter++;
         rsExtractor.extract( mol );
         if( rsExtractor.hasLargestRingSystem() )
         {  
            if( replaceOutputWithLargestRingSystem || replaceOutputWithLargestRing )
            {               
               // cache the current mol with sd tags
               OEGraphMol oldMol = new OEGraphMol();      
               oechem.OECopySDData(oldMol, mol);
               
               if ( replaceOutputWithLargestRingSystem )
               {
                  // Replace the current mol with the largest ring system mol
                  mol = rsExtractor.getLargestRingSystemMol();
               }
               else if ( replaceOutputWithLargestRing )
               {
                  // Replace the current mol with the largest ring 
                  mol = rsExtractor.getLargestRingMol();
               }
               
               // Save the SD data back...
               oechem.OECopySDData(mol, oldMol);
            }
            
            oechem.OESetSDData( mol, "largestRingSystem",
                                rsExtractor.getLargestRingSystemSMILES() );            
         }
         if( rsExtractor.getBasicRingSystemCount() > 0 )
         {   oechem.OESetSDData( mol, "basicRingSystems",
                                rsExtractor.getBasicRingSystemsSMILES() );
         }
         if( rsExtractor.hasLargestRing() )
         {   oechem.OESetSDData( mol, "largestRing",
                                rsExtractor.getLargestRingSMILES() );
         }
         oechem.OEWriteMolecule( outputOEThread, mol );

         //Output "." to show that the program is running.
         if( iCounter % 100 == 0 )
            System.err.print(".");
         if( iCounter % 4000 == 0 )
         {  System.err.printf( " %d %dsec\n",
                  iCounter, (System.currentTimeMillis()-start)/1000);
         }
         mol.Clear();
      }
      ifs.close();
      inFile = inFile.replaceAll( ".*" + Pattern.quote(File.separator), "" );
      System.err.printf( "SDFRingSystemExtractor: Read %d structures from %s. %d sec\n",
            iCounter, inFile, (System.currentTimeMillis()-start)/1000 );
   }


   private static void exitWithHelp(Options options) {
      HelpFormatter formatter = new HelpFormatter();
      formatter.printHelp( MY_NAME, options );
      System.exit(1);
   }


   /**
    * @param args
    */
   public static void main( String...args ) throws IOException
   {  // create command line Options object
      Options options = new Options();
      Option opt = new Option( OPT_INFILE, true, "input file oe-supported" );
      opt.setRequired( true );
      options.addOption( opt );

      opt = new Option( OPT_OUTFILE, true, "output file oe-supported" );
      opt.setRequired( false );
      options.addOption( opt );

      opt = new Option(OPT_REPLACE_OUTPUT_WITH_LARGEST_RING_SYSTEM, false, "If given, the largest ring system will be written to the output file instead of the input structure (default)");
      opt.setRequired(false);
      options.addOption(opt);
      
      opt = new Option(OPT_REPLACE_OUTPUT_WITH_LARGEST_RING, false, "If given, the largest single ring will be written to the output file instead of the input structure (default)");
      opt.setRequired(false);
      options.addOption(opt);
      
      opt = new Option(OPT_KEEP_ONE_SIDE_CHAIN_ATOM, false, "If given, the ring systems will retain one atom outside the ring.");
      opt.setRequired(false);
      options.addOption(opt);
      
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
      
      if( cmd.hasOption( OPT_REPLACE_OUTPUT_WITH_LARGEST_RING_SYSTEM ) )
      {  replaceOutputWithLargestRingSystem = true;}
         
      if( cmd.hasOption( OPT_REPLACE_OUTPUT_WITH_LARGEST_RING ) )
      {  replaceOutputWithLargestRing = true;}
      
      if( cmd.hasOption( OPT_KEEP_ONE_SIDE_CHAIN_ATOM ) )
      {  keepOneSideChainAtom = true;}
      
      SDFRingSystemExtractor extractor = new SDFRingSystemExtractor( outFile );
      extractor.run( inFile );
      extractor.close();
   }
}
