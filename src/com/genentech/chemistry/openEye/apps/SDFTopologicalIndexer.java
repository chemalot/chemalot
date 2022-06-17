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
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.TimeoutException;
import java.util.regex.Pattern;

import openeye.oechem.*;

import org.apache.commons.cli.*;

import com.aestel.utility.DataFormat;
import com.genentech.chemistry.openEye.topoIndexes.TopologicalIndexer;
import com.genentech.oechem.tools.OETools;




/**

 * @author Alberto Gobbi/ 2012
 * Copyright 2012 Genentech
 */
public class SDFTopologicalIndexer
{  private static final String MY_NAME = "SDFTopologicalIndexer";

   private static final String OPT_INFILE    = "in";
   private static final String OPT_OUTFILE   = "out";

   private final oemolothread outputOEThread;
   private final TopologicalIndexer indexer;

   private static final Set<String> AVAILIndexes = new HashSet<String>();

   private final boolean doJ;
   private final boolean doJStar;
   private final boolean doJX;
   private final boolean doJXStar;
   private final boolean doJY;
   private final boolean doJYStar;
   private final boolean doWiener;
   private final boolean doZagreb;


   static
   {  AVAILIndexes.add("J");
      AVAILIndexes.add("JStar");
      AVAILIndexes.add("JX");
      AVAILIndexes.add("JXStar");
      AVAILIndexes.add("JY");
      AVAILIndexes.add("JYStar");
      AVAILIndexes.add("Wiener");
      AVAILIndexes.add("Zagreb");
   }

   private SDFTopologicalIndexer( String outFile, Set<String> selectedIndexes )
   {  outputOEThread = new oemolothread(outFile);
      indexer        = new TopologicalIndexer(  );

      doJ            = selectedIndexes.contains("J");
      doJStar        = selectedIndexes.contains("JStar");
      doJX           = selectedIndexes.contains("JX");
      doJXStar       = selectedIndexes.contains("JXStar");
      doJY           = selectedIndexes.contains("JY");
      doJYStar       = selectedIndexes.contains("JYStar");
      doWiener       = selectedIndexes.contains("Wiener");
      doZagreb       = selectedIndexes.contains("Zagreb");
   }


   private void close()
   {  outputOEThread.close();
      indexer.close();
   }



   private void run( String inFile )
   {  oemolithread ifs = new oemolithread(inFile);
      long start = System.currentTimeMillis();
      int iCounter = 0; //Structures in the SD file.

      OEMolBase mol = new OEGraphMol();
      while( oechem.OEReadMolecule( ifs, mol ) )
      {  iCounter++;

         try
         {  indexer.computeIndexes(mol);
//System.err.println(OETools.molToCanSmi(mol, true));
            if( doJ )
               oechem.OESetSDData(mol, "J",     DataFormat.formatNumber(indexer.getBalabanJIndex(), "r3"));
            if( doJStar )
               oechem.OESetSDData(mol, "JStar", DataFormat.formatNumber(indexer.getBalabanJStarIndex(), "r3"));
            if( doJX )
               oechem.OESetSDData(mol, "JX",    DataFormat.formatNumber(indexer.getBalabanJXIndex(), "r3"));
            if( doJXStar )
               oechem.OESetSDData(mol, "JXStar",DataFormat.formatNumber(indexer.getBalabanJXStarIndex(), "r3"));
            if( doJY )
               oechem.OESetSDData(mol, "JY",    DataFormat.formatNumber(indexer.getBalabanJYIndex(), "r3"));
            if( doJYStar )
               oechem.OESetSDData(mol, "JYStar",DataFormat.formatNumber(indexer.getBalabanJYStarIndex(), "r3"));
            if( doWiener )
               oechem.OESetSDData(mol, "Wiener",Long.toString(indexer.getWienerIndex()));
            if( doZagreb )
               oechem.OESetSDData(mol, "Zagreb",Integer.toString(indexer.getZagrebIndex()));

         } catch(TimeoutException e)
         {  System.err.println("Conmputing topological Index timed out for: " + OETools.molToCanSmi(mol, true));
         } catch(NumberFormatException e)
         {  System.err.println("Conmputing topological Index caused NumberFormatException for: " + OETools.molToCanSmi(mol, true));
            e.printStackTrace(System.err);
         } catch(OutOfMemoryError e)
         {  System.err.println("Conmputing topological Index caused OutOfMemoryError for: " + OETools.molToCanSmi(mol, true));
            throw e;
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
      ifs.close();
      inFile = inFile.replaceAll( ".*" + Pattern.quote(File.separator), "" );
      System.err.printf( "%s: Read %d structures from %s. in %d sec\n",
               MY_NAME, iCounter, inFile, (System.currentTimeMillis()-start)/1000 );
   }





   private static void exitWithHelp(Options options)
   {  int width = 100;

      StringBuilder types = new StringBuilder();
      types.append("'all' or one or multiple of: ");
      for( String typ : AVAILIndexes ) types.append(typ).append(" ");

      String start = MY_NAME + " all|IndexTypes";
      String head  = "Will generate topological indexes for input compounds.";

      HelpFormatter formatter = new HelpFormatter();
      formatter.setArgName("p");
      PrintWriter out = new PrintWriter(System.err,true);
      formatter.printUsage(out, width, start, options);
      formatter.printWrapped(out, width, head);
      formatter.printOptions(out, width, options, 2, 2);
      formatter.printWrapped(out, width, 18, "  IndexTypes:     " + types.toString());
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
      opt.setArgName("fn");
      options.addOption( opt );

      opt = new Option( OPT_OUTFILE, true,
               "output file oe-supported. Use .sdf|.smi to specify the file type." );
      opt.setRequired( true );
      opt.setArgName("fn");
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

      if( args.length == 0 )
      {  System.err.println("Specify at least one index type");
         exitWithHelp(options);
      }

      String inFile = cmd.getOptionValue( OPT_INFILE );
      String outFile = cmd.getOptionValue( OPT_OUTFILE );
      Set<String> selectedIndexes = new HashSet<String>(args.length);
      if( args.length == 1 && "all".equalsIgnoreCase(args[0]))
         selectedIndexes = AVAILIndexes;
      else
         selectedIndexes.addAll(Arrays.asList(args));

      if( ! AVAILIndexes.containsAll(selectedIndexes) )
      {  selectedIndexes.removeAll(AVAILIndexes);
         StringBuilder err = new StringBuilder("Unknown Index types: ");
         for( String it : selectedIndexes ) err.append(it).append(" ");
         System.err.println(err);
         exitWithHelp(options);
      }

      SDFTopologicalIndexer sdfIndexer = new SDFTopologicalIndexer(outFile, selectedIndexes);

      sdfIndexer.run( inFile );
      sdfIndexer.close();
   }
}
