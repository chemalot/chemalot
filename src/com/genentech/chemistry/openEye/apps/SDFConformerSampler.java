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

import com.genentech.chemistry.openEye.conformerSampler.ConformerSampler;




/**

 * @author Alberto Gobbi/ 2012
 * Copyright 2012 Genentech
 */
public class SDFConformerSampler
{  private static final String MY_NAME = "SDFConformerSampler";

   private static final String OPT_INFILE    = "in";
   private static final String OPT_OUTFILE   = "out";
   private static final String OPT_TORSION_FILE = "torsionFile";
   private static final String OPT_MAX_CONFS    = "maxConfs";

   private final oemolothread outputOEThread;
   private final ConformerSampler scanner;

   private final long maxConfs;

   private SDFConformerSampler( String smartsFile, String outFile, long maxConfs )
   {  this.maxConfs  = maxConfs;
      outputOEThread = new oemolothread(outFile);
      scanner        = new ConformerSampler( smartsFile );
   }


   private void close()
   {  outputOEThread.close();
      scanner.close();
   }



   private void run( String inFile )
   {  oemolithread ifs = new oemolithread(inFile);
      long start = System.currentTimeMillis();
      int iCounter = 0; //Structures in the SD file.
      int oCounter = 0;

      OEMolBase mol = new OEGraphMol();
      while( oechem.OEReadMolecule( ifs, mol ) )
      {  iCounter++;

         OEMCMolBase confs = scanner.createConformations(mol, maxConfs);
         oechem.OEWriteMolecule(outputOEThread, confs);
         oCounter += confs.NumConfs();
         confs.delete();

         //Output "." to show that the program is running.
         if( iCounter % 100 == 0 )
            System.err.print(".");
         if( iCounter % 4000 == 0 )
         {  System.err.printf( " %d %dsec\n",
                  iCounter, (System.currentTimeMillis()-start)/1000);
         }
      }

      mol.delete();
      ifs.close();
      inFile = inFile.replaceAll( ".*" + Pattern.quote(File.separator), "" );
      System.err.printf( "%s: Read %d structures from %s. Written %d confs in %d sec\n",
               MY_NAME, iCounter, inFile, oCounter, (System.currentTimeMillis()-start)/1000 );
   }





   private static void exitWithHelp(Options options)
   {  HelpFormatter formatter = new HelpFormatter();
      String head = "Will generate maxConf conformers by rotating torsional angels as "
                   +"defined in the torsion file. The default torsion file will rotate "
                   +"OH and NH2 groups. The input conformation is always returned first.";
      formatter.printHelp( MY_NAME, head, options, "", true );
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

      opt = new Option( OPT_MAX_CONFS, true,
               "Maximum number of conformations per input." );
      opt.setRequired( true );
      options.addOption( opt );

      opt = new Option( OPT_TORSION_FILE, true,
               "Optional: to overwrite torsion definition file." );
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
      String smartsFile = cmd.getOptionValue( OPT_TORSION_FILE );
      long maxConfs = Long.parseLong(cmd.getOptionValue( OPT_MAX_CONFS ));

      SDFConformerSampler scanner = new SDFConformerSampler(smartsFile, outFile, maxConfs);

      scanner.run( inFile );
      scanner.close();
   }
}
