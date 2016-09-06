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
import java.util.EnumSet;
import java.util.regex.Pattern;

import openeye.oechem.*;

import org.apache.commons.cli.*;

import com.genentech.chemistry.openEye.cats.AtomTyperInterface;
import com.genentech.chemistry.openEye.cats.CATSIndexer;
import com.genentech.chemistry.openEye.cats.CATSIndexer.Normalization;




/**

 * @author Alberto Gobbi/ 2012
 * Copyright 2012 Genentech
 */
public class SDFCatsIndexer
{  private static final String MY_NAME = "SDFCatsIndexer";

   private static final String OPT_INFILE    = "in";
   private static final String OPT_OUTFILE   = "out";
   private static final String OPT_NORMALIZATION = "normalization";
   private static final String OPT_PRINTDESC = "printDescriptors";
   private static final String OPT_RGROUPTYPES = "rGroups";

   private final CATSIndexer indexer;

   private SDFCatsIndexer(AtomTyperInterface[] myTypes, String tagPrefix)
   {  indexer = new CATSIndexer(myTypes, tagPrefix);
   }


   private void close()
   {  indexer.close();
   }



   private void run( String inFile, String outFile, EnumSet<CATSIndexer.Normalization> normMeth )
   {  oemolithread ifs = new oemolithread(inFile);
      oemolothread ofs = new oemolothread(outFile);

      long start = System.currentTimeMillis();
      int iCounter = 0; //Structures in the SD file.

      OEMolBase mol = new OEGraphMol();
      while( oechem.OEReadMolecule( ifs, mol ) )
      {  iCounter++;

         indexer.compute2DCats(mol, normMeth);

         oechem.OEWriteMolecule(ofs, mol);

         //Output "." to show that the program is running.
         if( iCounter % 100 == 0 )
            System.err.print(".");
         if( iCounter % 4000 == 0 )
         {  System.err.printf( " %d %dsec\n",
                  iCounter, (System.currentTimeMillis()-start)/1000);
         }
      }

      mol.delete();
      ofs.close();
      ifs.close();
      ofs.delete();
      ifs.delete();
      inFile = inFile.replaceAll( ".*" + Pattern.quote(File.separator), "" );
      System.err.printf( "%s: Read %d structures from %s. in %d sec\n",
               MY_NAME, iCounter, inFile, (System.currentTimeMillis()-start)/1000 );
   }

   private void printDescriptors( String inFile, String outFile ) throws FileNotFoundException
   {  oemolithread ifs = new oemolithread(inFile);
      PrintStream out;
      if( ".tab".equals(outFile))
         out = System.out;
      else
         out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outFile)));

      long start = System.currentTimeMillis();
      int iCounter = 0; //Structures in the SD file.

      indexer.printDescriptorHeader(out);
      OEMolBase mol = new OEGraphMol();
      while( oechem.OEReadMolecule( ifs, mol ) )
      {  iCounter++;

         indexer.printDistanceDescriptors(out, mol);

         //Output "." to show that the program is running.
         if( iCounter % 100 == 0 )
            System.err.print(".");
         if( iCounter % 4000 == 0 )
         {  System.err.printf( " %d %dsec\n",
                  iCounter, (System.currentTimeMillis()-start)/1000);
         }
      }

      mol.delete();
      out.close();
      ifs.close();
      ifs.delete();
      inFile = inFile.replaceAll( ".*" + Pattern.quote(File.separator), "" );
      System.err.printf( "%s: Read %d structures from %s. in %d sec\n",
               MY_NAME, iCounter, inFile, (System.currentTimeMillis()-start)/1000 );
   }





   private static void exitWithHelp(Options options)
   {  int width = 100;

      String start = MY_NAME;
      String head  = "Will generate CATS Fingerprints for input compounds.";

      HelpFormatter formatter = new HelpFormatter();
      formatter.setArgName("p");
      PrintWriter out = new PrintWriter(System.err,true);
      formatter.printUsage(out, width, start, options);
      formatter.printWrapped(out, width, head);
      formatter.printOptions(out, width, options, 2, 2);
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

      opt = new Option( OPT_NORMALIZATION, true,
               "Normalization method: Counts|CountsPerAtom|CountsPerFeature(def) multiple allowed" );
      opt.setArgName("meth");
      options.addOption( opt );

      opt = new Option( OPT_PRINTDESC, false,
               "Causes the descriptor for describing each linear path in a molceule to be created");
      options.addOption( opt );

      opt = new Option( OPT_RGROUPTYPES, false,
               "treat RGroup attachement point ([U]) as atom type.");
      options.addOption( opt );

      CommandLineParser parser = new PosixParser();
      CommandLine cmd = null;
      try
      {  cmd = parser.parse( options, args );
      } catch( Exception e )
      {  System.err.println( e.getMessage() );
         exitWithHelp( options );
      }

      String inFile = cmd.getOptionValue( OPT_INFILE );
      String outFile = cmd.getOptionValue( OPT_OUTFILE );

      AtomTyperInterface[] myTypes = CATSIndexer.typers;
      String tagPrefix = "";
      if( cmd.hasOption(OPT_RGROUPTYPES))
      {  myTypes = CATSIndexer.rgroupTypers;
         tagPrefix = "RG";
      }

      if( cmd.hasOption(OPT_PRINTDESC) )
      {  SDFCatsIndexer sdfIndexer = new SDFCatsIndexer(myTypes, tagPrefix);
         sdfIndexer.printDescriptors(inFile, outFile);
         sdfIndexer.close();
         return;
      }

      EnumSet<Normalization> normMeth = EnumSet.noneOf(Normalization.class);
      if( cmd.hasOption(OPT_NORMALIZATION) )
         for(String n: cmd.getOptionValues(OPT_NORMALIZATION))
            normMeth.add(Normalization.valueOf(n));
      else
         normMeth.add(Normalization.CountsPerFeature);

      SDFCatsIndexer sdfIndexer = new SDFCatsIndexer(myTypes, tagPrefix);
      sdfIndexer.run( inFile, outFile, normMeth );
      sdfIndexer.close();
   }
}
