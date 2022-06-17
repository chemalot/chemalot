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

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import openeye.oechem.*;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.PosixParser;

import com.aestel.utility.DataFormat;

/**
 * Calculate the RMSD similarity of two molecule (docking) poses
 *
 * Output will also include avg(RMSD) and/or median(RMSD) on per compound basis
 *
 * @author Man-Ling Lee / March 04, 2015
 * Copyright 2015 Genentech, Inc.
 */
public class SdfRMSDNNFinder
{  private static final String MY_NAME = "sdfRMSDNNFinder";

   private static final String OPT_INFILE   = "in";
   private static final String OPT_OUTFILE  = "out";
   private static final String OPT_REFFILE  = "ref";
   private static final String OPT_MOLIdTag = "molIdTag";
   private static final String OPT_MIRROR   = "mirror";
   private static final String OPT_DONotOpt = "doNotOptimize";

   private static final String OTAG_MOL1 = "Pose_Idx";
   private static final String OTAG_MOL2 = "NN_Pose_Idx";
   private static final String OTAT_RMSD = "NN_RMSD";

   private final String  molIdTag;
   private final boolean doMirror;
   private final boolean doOptimize;

   private final Map<String,List<OEMolBase>> refMolPoses;
   private final oemolothread outputOEThread;
   private final boolean noRefFile;


   private SdfRMSDNNFinder( String refFile, String outFile, String molIdTag,
               boolean doMirror, boolean doOptimize, boolean noRefFile )
   {  outputOEThread = new oemolothread( outFile );
      this.molIdTag = molIdTag;
      this.doMirror = doMirror;
      this.doOptimize = doOptimize;
      this.noRefFile = noRefFile;

      refMolPoses = new HashMap<String,List<OEMolBase>>();
      OEMolBase mol = new OEGraphMol();
      oemolithread refOEThread = new oemolithread( refFile );
      while( oechem.OEReadMolecule( refOEThread, mol ) )
      {
         String molId = oechem.OEGetSDData( mol, molIdTag );
         List<OEMolBase> list = refMolPoses.get( molId );
         if( list == null )
         { list = new ArrayList<OEMolBase>();
            refMolPoses.put( molId, list );
         }
         oechem.OEAssignAromaticFlags(mol);
         list.add( new OEGraphMol(mol) );
      }
      mol.delete();
      refOEThread.close();
      refOEThread.delete();
   }


   private void run( String inFile )
   {  oemolithread inputOEThread = new oemolithread( inFile );
      long start = System.currentTimeMillis();
      int iCounter = 0; //records in the input file.
      int oCounter = 0; //records in the output file, i.e. similarity pair

      OEMolBase inPose  = new OEGraphMol();
      int myPoseIndex = -1;
      String previousMolId = "";
      while( oechem.OEReadMolecule( inputOEThread, inPose ) )
      {
         iCounter++;
         //Output "." to show that the program is progressing
         if( iCounter % 100 == 0 )
            System.err.print( "." );
         if( iCounter % 4000 == 0 )
            System.err.printf( " %d %dsec\n",
                     iCounter, (System.currentTimeMillis()-start)/1000 );

         String molId = oechem.OEGetSDData( inPose, molIdTag );
         if( molId.equals( previousMolId ) )
         {  ++myPoseIndex;
         } else
         {  myPoseIndex = 0;
            previousMolId = molId;
         }
         oechem.OESetSDData( inPose, OTAG_MOL1, Integer.toString( myPoseIndex ) );

         List<OEMolBase> refPoses = refMolPoses.get( molId );
         if( refPoses == null )
         {  oechem.OEWriteMolecule( outputOEThread, inPose );
            continue;
         }

         // apply same aromaticity problems as on reference mols
         oechem.OEAssignAromaticFlags(inPose);

         OEMolBase mirrorPose = null;
         if( doMirror && ! SdfRMSDSphereExclusion.isChiral( inPose ) )
         {  mirrorPose = new OEGraphMol( inPose );
            SdfRMSDSphereExclusion.createMirror( mirrorPose );
         }

         for( int i=0; i<refPoses.size(); i++ )
         {
            /*Do not compare a pose with itself in case inFile is used as refFile
             *as well. If refFile is specified, it is assumed that refFile contains
             *poses from different docking run. Therefore, every poses in inFile
             *are compared with every poses in the refFile having the same molecule
             *identifier.*/
            if( noRefFile && i == myPoseIndex )
            {  if( refPoses.size() == 1 )
                  oechem.OEWriteMolecule( outputOEThread, inPose );
               continue;
            }
            OEMolBase refPose = refPoses.get( i );

            double rmsd = oechem.OERMSD( inPose, refPose, false, true, doOptimize );
            if( rmsd < 0D )
            {  System.err.println( "OERMSD returned -1. " +
                        "Are you comparing two different structures?" );
               continue;
            }
            if( mirrorPose != null )
            {  double mirrorRmsd = oechem.OERMSD(
                        mirrorPose, refPose, true, true, doOptimize );
               if( mirrorRmsd < rmsd )
                  rmsd = mirrorRmsd;
            }
            // Writing inPose directly with OEWriteMolecule causes OEBug with OERMSD
            // returning -1.0 for all comparisons following the first
            oechem.OESetSDData( inPose, OTAG_MOL2, Integer.toString( i ) );
            oechem.OESetSDData( inPose, OTAT_RMSD, DataFormat.formatNumber(rmsd, "si3") );
            oechem.OEWriteMolecule( outputOEThread, inPose );
            ++oCounter;
         }
         if( mirrorPose != null )
            mirrorPose.delete();
      }
      inPose.delete();
      inputOEThread.close();
      inputOEThread.delete();
      System.err.printf( "\n%s: Read %d molecules. Output %d molecule pairs in %d sec\n\n",
               MY_NAME, iCounter, oCounter, (System.currentTimeMillis()-start)/1000 );
   }


   private static String writeRefMolPosesToTempFile( String fileName )
   throws IOException
   {  String tmpFileName = File.createTempFile( "rmsdNN", ".oeb" )
                               .getAbsolutePath();
      OEMolBase mol    = new OEGraphMol();
      oemolothread out = new oemolothread( tmpFileName );
      oemolithread in  = new oemolithread( fileName );

      while( oechem.OEReadMolecule( in, mol ) )
         oechem.OEWriteMolecule( out, mol );

      mol.delete();
      out.close();
      out.delete();
      in.close();
      in.delete();
      return tmpFileName;
   }


   private void close()
   {  outputOEThread.close();
      outputOEThread.delete();

      Iterator<String> keys = refMolPoses.keySet().iterator();
      while( keys.hasNext() )
      {  List<OEMolBase> poses = refMolPoses.get( keys.next() );
         for( OEMolBase pose : poses )
            pose.delete();
      }
      refMolPoses.clear();
   }


   private static void exitWithHelp( Options options )
   {  HelpFormatter formatter = new HelpFormatter();
      String head = "Calculate the pairwise RMSD similarity of two docking poses.\n"
                  + "The program compares records where the molID of the reference matches the input.\n"
                  + "If reference file is specified, program compares the poses in\n"
                  + "the input file against poses of the same molecule in the\n"
                  + "reference file. Otherwise the input compounds are used as reference.\n"
                  + "The output will include the following tags:\n"
                  + "   Pose_Idx the index of the input record with same molID in order of RMSD\n"
                  + "   NN_Pose_Idx the index of the reference recordwith same molID matching\n"
                  + "   NN_RMSD the RMSD\n";
      formatter.printHelp( MY_NAME, head, options, "", true );
      System.exit(1);
   }


   public static void main( String[] args ) throws IOException
   {
      Options options = new Options();
      Option opt = new Option( OPT_INFILE, true,
               "input file oe-supported Use .sdf|.smi to specify the file type." );
      opt.setRequired( true );
      options.addOption( opt );

      opt = new Option( OPT_OUTFILE, true,
               "output file oe-supported. Use .sdf|.smi to specify the file type." );
      opt.setRequired( true );
      options.addOption( opt );

      opt = new Option( OPT_REFFILE, true,
               "Reference file containing poses from reference docking run. " +
               "If " + OPT_REFFILE + " not specified, program uses the input file " +
               "(internal NN analysis)" );
      options.addOption( opt );

      opt = new Option( OPT_MOLIdTag, true,
               "Name of the field containing the molecule identifier. " +
               "Only conformations with same ID are compared." +
               " Assumption: Ref file uses the same field name.");
      opt.setRequired( true );
      options.addOption( opt );

      opt = new Option( OPT_DONotOpt, false,
               "If specified the RMSD is computed without trying to optimize the alignment." );
      options.addOption( opt );

      opt = new Option( OPT_MIRROR,false,
               "For non-chiral molecules also try mirror image" );
      options.addOption( opt );

      CommandLineParser parser = new PosixParser();
      CommandLine cmd = null;
      try
      {  cmd = parser.parse( options, args );
      } catch ( Exception e )
      {  System.err.println( e.getMessage() );
         exitWithHelp( options );
      }
      args = cmd.getArgs();
      if( args.length != 0 )
      {
         ; //error message
      }

      String inFile = cmd.getOptionValue( OPT_INFILE );
      String outFile = cmd.getOptionValue( OPT_OUTFILE );

      String refFile = cmd.getOptionValue( OPT_REFFILE );
      boolean noRefFile = false;
      if( refFile == null )
      {  if( inFile.startsWith( "." ) )
            inFile = writeRefMolPosesToTempFile( inFile );
         refFile = inFile;
         noRefFile = true;
      }

      String molIdTag = cmd.getOptionValue( OPT_MOLIdTag );
      boolean doOptimize = ! cmd.hasOption( OPT_DONotOpt );
      boolean doMirror = cmd.hasOption( OPT_MIRROR );

      SdfRMSDNNFinder nnFinder = new SdfRMSDNNFinder(
               refFile, outFile, molIdTag, doMirror, doOptimize, noRefFile );

      nnFinder.run( inFile );
      nnFinder.close();
   }

}
