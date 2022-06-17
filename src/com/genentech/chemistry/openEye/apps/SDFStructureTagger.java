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
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map.Entry;
import java.util.regex.Pattern;

import openeye.oechem.OEGraphMol;
import openeye.oechem.OEIFlavor;
import openeye.oechem.OEMDLQueryOpts;
import openeye.oechem.OEMolBase;
import openeye.oechem.OEQMol;
import openeye.oechem.OESubSearch;
import openeye.oechem.oechem;
import openeye.oechem.oemolistream;
import openeye.oechem.oemolithread;
import openeye.oechem.oemolothread;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.PosixParser;

import com.aestel.utility.Message;
import com.aestel.utility.Message.Level;
import com.aestel.utility.exception.IncorrectInputException;
import com.genentech.chemistry.openEye.StructureTagger;




/**
 * Use it to tag a structure based on SMARTS or molfile (queries) input.
 *
 * @author Man-Ling Lee / February 12, 2011
 * Copyright 2011-2012 Genentech
 */
public class SDFStructureTagger
{  private static final String MY_NAME = "SDFStructureTagger";
   private static final String OPT_INFILE    = "in";
   private static final String OPT_OUTFILE   = "out";
   private static final String OPT_SMARTS_FILE = "smarts";
   private static final String OPT_QSDF_FILE = "sdfile";
   private static final String OPT_TAG_SETS = "sets";
   private static final String OPT_TAG_INFO = "tag_info";
   private static final String OPT_OUTPUT_COUNT = "output_count";
   private static final String OPT_OUTPUT_EXISTS = "output_exists";
   private static final String OPT_ONLY_MATCHES  = "onlyMatches";
   private static final String OPT_NO_MATCHES  = "noMatches";

   public static enum InputFieldName
   {  TAG_NAME,
      SET_NAME;
      public static InputFieldName toInputFieldName( String str )
      {  assert str != null : "Must specify str.";
         try
         {  return InputFieldName.valueOf( str.toUpperCase() );
         } catch( IllegalArgumentException iae )
         { return null; }
      }
   }

   private static final String TAG_SETS_ALL = "all";

   public static enum OutputTag
   {  tagCount,
      allTags,
      allTagPattern,
      firstTag,
      firstTagPattern;
      public static OutputTag toOutputTag( String str )
      {  OutputTag[] outTags = OutputTag.values();
         for( int i=0; i<outTags.length; i++ )
         {  if( outTags[i].toString().equalsIgnoreCase( str ) )
               return outTags[i];
         }
         return null;
      }
      public static String getTagsAsText( String separator )
      {  assert separator != null && separator.length() > 0 :
               "Must specify a separator.";
         OutputTag[] tags = OutputTag.values();
         StringBuffer sb = new StringBuffer().append( tags[0] );
         for( int i=1; i<tags.length; i++ )
            sb.append( separator ).append( tags[i] );
         return sb.toString();
      }
   }


   private final oemolothread outputOEThread;
   private StructureTagger structureTagger;
   private HashSet<OutputTag> requestedOutFields;
   private boolean outputOccurrence;
   private boolean outputExists;
   private boolean onlyMatches;
   private boolean noMatches;
   //private List<String> requestedTagFields;


   private SDFStructureTagger( String outFile )
   {  outputOEThread     = new oemolothread(outFile);
      structureTagger    = new StructureTagger();
      requestedOutFields = new HashSet<OutputTag>();
   }


   private void close()
   {  outputOEThread.close();
   }


   private void prepare( String smartsFile, String sdfFile, String setsForTagging,
            String requestedTagInfo, boolean outputOccurrence, boolean outputExists, boolean onlyMatches, boolean noMatches )
   throws IncorrectInputException, IOException
   {  if( ( smartsFile == null || smartsFile.length() == 0 ) &&
          ( sdfFile == null || sdfFile.length() == 0 ) )
      {  String s = "Both " + OPT_SMARTS_FILE + " and " + OPT_QSDF_FILE + " are not set.";
         throw new IncorrectInputException( new Message( s, Level.ERROR, null ) );
      }

      //Prepare the structure matching pattern for tagging
      String[] setNames = setsForTagging.split( "\\|" );
      for( int i=0; i<setNames.length; i++ )
         setNames[i] = setNames[i].trim();

      List<Message> messageList = new ArrayList<Message>();
      try
      {  readSMARTS( smartsFile, setNames );
      } catch( FileNotFoundException e )
      {  String s = "Cannot open " + smartsFile + ".\n";
         messageList.add( new Message( s, Level.ERROR, null ) );
      }
      try
      {  readSDFile( sdfFile, setNames );
      } catch( FileNotFoundException e )
      {  String s = "Cannot open " + sdfFile + ".\n";
         messageList.add( new Message( s, Level.ERROR, null ) );
      }

      //Validate OPT_OUT_TAGS
      if( requestedTagInfo != null && requestedTagInfo.length() > 0 )
      {  String[] tags = requestedTagInfo.split( "\\|" );
         for( int i=0; i<tags.length; i++ )
         {  String tag = tags[i].trim();
            if( OutputTag.toOutputTag( tag ) != null )
            {  requestedOutFields.add( OutputTag.toOutputTag( tag ) );
            } else
            {  String s = tag + ": Incorrect " + OPT_TAG_INFO + " flag.\n";
               messageList.add( new Message( s, Level.ERROR, null ) );
            }
         }
      }
      this.outputOccurrence = outputOccurrence;
      this.outputExists = outputExists;
      this.onlyMatches = onlyMatches;
      this.noMatches = noMatches;
   }


   /**
    * Parse the tab-delimited file with the SMARTS definitions and add them to
    * the structureTagger class instance.
    *
    * @param file Name should include the path.
    * @param keepSetNames  specify with sets to be used for tagging. If keepSetNames
    *          is null or does not contain any element, all sets will be used.
    *
    * @return  the number of smarts retained for tagging.
    *
    * @throws FileNotFoundException
    * @throws IOException
    */
   private int readSMARTS( String file, String[] keepSetNames )
   throws FileNotFoundException, IOException
   {  if( file == null || file.length() == 0 )
         return 0;

      BufferedReader reader = new BufferedReader( new FileReader( file ) );
      String line;
      int smartsCounter = 0;

      while( null != ( line = reader.readLine() ) )
      {  if( line.length() == 0 ||
             line.startsWith( "#" ) || line.startsWith( "\"#" ) )
            continue;

         List<String> falseLines = new ArrayList<String>();
         String[] fields = line.split( "\t" ); //SMARTS, tag name, set name
         if( fields.length < 3 )
         {  falseLines.add( line );
            continue;
         }
         String smarts = fields[0].trim().replaceAll( "\"", "" );
         String tag    = fields[1].trim().replaceAll( "\"", "" );
         String set    = fields[2].trim().replaceAll( "\"", "" );
         if( keepSet( fields[2].trim(), keepSetNames ) )
         {  OESubSearch subSearch = new OESubSearch( smarts );
            structureTagger.addPattern( smarts, subSearch, tag, set );
            smartsCounter++;
         }
      }
      reader.close();
      System.err.printf( "Read %d smarts\n", smartsCounter );
      return smartsCounter;
   }


   private static boolean keepSet( String inSetName, String[] setsForTagging )
   {  if( setsForTagging[0].contains( TAG_SETS_ALL ) )
         return true;

      for( int i=0; i<setsForTagging.length; i++ )
      {  if( setsForTagging[i].equalsIgnoreCase( inSetName ) )
            return true;
      }
      return false;
   }


   /**
    * Parse the SD file with the MDL query definitions and add them to the
    * structureTagger class instance.
    *
    * @param file Name should include the path.
    * @param keepSetNames  specify with sets to be used for tagging. If keepSetNames
    *          is null or does not contain any element, all sets will be used.
    *
    * @return  the number of MDL query definitions retained for tagging.
    *
    * @throws FileNotFoundException
    * @throws IOException
    */
   private int readSDFile( String file, String[] keepSetNames )
   throws FileNotFoundException, IOException
   {  if( file == null || file.length() == 0 )
         return 0;

      // read queries into List
      oemolistream qfile = new oemolistream( file );
      int aromodel = OEIFlavor.Generic.OEAroModelOpenEye;
      int qflavor  = qfile.GetFlavor( qfile.GetFormat() );
      qfile.SetFlavor( qfile.GetFormat(), ( qflavor|aromodel ) );
      int opts = OEMDLQueryOpts.Default|OEMDLQueryOpts.SuppressExplicitH;
      OEGraphMol mol = new OEGraphMol();
      OEQMol qmol = new OEQMol();
      int nMol = 0;
      
      while( oechem.OEReadMDLQueryFile( qfile, mol ) )
      {  String tagName = oechem.OEGetSDData( mol,
                  InputFieldName.TAG_NAME.toString() ).trim();
         String setName = oechem.OEGetSDData( mol,
                  InputFieldName.SET_NAME.toString() ).trim();
         
         oechem.OEBuildMDLQueryExpressions(qmol,mol, opts);
         OESubSearch subSearch = new OESubSearch( qmol );
         
         if( keepSet( setName, keepSetNames ) )
         {  structureTagger.addPattern( null, subSearch, tagName, setName );
            nMol++;
         }
         mol.Clear();
         qmol.Clear();
      }
      System.err.printf("Read %d query molecules\n", nMol);

      mol.delete();
      qmol.delete();
      qfile.close();
      qfile.delete();
      return nMol;
   }


   private void run( String inFile )
   {  oemolithread ifs = new oemolithread( inFile );
      long start = System.currentTimeMillis();
      int iCounter = 0; //Structures in the SD file.

      OEMolBase mol = new OEGraphMol();
      while( oechem.OEReadMolecule( ifs, mol ) )
      {  iCounter++;

         boolean isMatch = false;
         if( ( requestedOutFields.size() == 0 && !outputOccurrence && !outputExists )
          || requestedOutFields.size() > 0 )
            isMatch |= addTaggingInfoToOutput( mol );
         if( outputOccurrence )
            isMatch |= addOccurrenceToOutput( mol );
         if( outputExists )
            isMatch |= addExistsToOutput( mol );

         if( (! onlyMatches || isMatch) && !noMatches )
            oechem.OEWriteMolecule( outputOEThread, mol );

         if( noMatches && !isMatch && !onlyMatches )
            oechem.OEWriteMolecule( outputOEThread, mol );
         
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
      ifs.delete();
      inFile = inFile.replaceAll( ".*" + Pattern.quote(File.separator), "" );
      System.err.printf( "%s: Read %d structures from %s. %d sec\n",
               MY_NAME, iCounter, inFile, (System.currentTimeMillis()-start)/1000 );
   }


   /**
    *
    * @return true if at least one match occured.
    */
   private boolean addOccurrenceToOutput( OEMolBase mol )
   {  LinkedHashMap<String,Integer> map = structureTagger.countOccurrence( mol );
      Iterator<Entry<String,Integer>> iterator = map.entrySet().iterator();
      while( iterator.hasNext() )
      {  Entry<String,Integer> entry = iterator.next();
         oechem.OESetSDData( mol, entry.getKey(), String.valueOf( entry.getValue() ) );
      }
      return map.size() > 0;
   }

   /**
    *
    * @return true if at least one match occured.
    */
   private boolean addExistsToOutput( OEMolBase mol )
   {  LinkedHashMap<String,Boolean> map = structureTagger.checkOccurrence( mol );
      Iterator<Entry<String,Boolean>> iterator = map.entrySet().iterator();
      while( iterator.hasNext() )
      {  Entry<String,Boolean> entry = iterator.next();
         oechem.OESetSDData( mol, entry.getKey(),
                  entry.getValue().booleanValue() ? "1" : "0" );
      }
      return map.size() > 0;
   }


   private boolean addTaggingInfoToOutput( OEMolBase mol )
   {  boolean matched = structureTagger.tagStructure( mol );

      if( requestedOutFields.size() == 0 ||
          requestedOutFields.contains( OutputTag.tagCount ) )
      {  oechem.OESetSDData( mol, OutputTag.tagCount.toString(),
                  String.valueOf( structureTagger.getTagCounts() ) );
      }
      if( requestedOutFields.size() == 0 ||
               requestedOutFields.contains( OutputTag.allTags ) )
      {  oechem.OESetSDData( mol, OutputTag.allTags.toString(),
                  structureTagger.getAllTags() );
      }
      if( requestedOutFields.size() == 0 ||
               requestedOutFields.contains( OutputTag.allTagPattern ) )
      {  oechem.OESetSDData( mol, OutputTag.allTagPattern.toString(),
                  structureTagger.getAllTagPattern() );
      }
      if( requestedOutFields.size() == 0 ||
               requestedOutFields.contains( OutputTag.firstTag ) )
      {  oechem.OESetSDData( mol, OutputTag.firstTag.toString(),
                  structureTagger.getFirstTag() );
      }
      if( requestedOutFields.size() == 0 ||
               requestedOutFields.contains( OutputTag.firstTagPattern ) )
      {  oechem.OESetSDData( mol, OutputTag.firstTagPattern.toString(),
                  structureTagger.getFirstTagPattern() );
      }

      return matched;
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
               "Tag-delimited file containing SMARTS, tag name, name of the tag set." +
               " Lines starting with '#' are ignores as comment." +
               " NOTE: Use the SMARTS features supported by OpenEye." );
      opt.setRequired( false );
      options.addOption( opt );

      opt = new Option( OPT_QSDF_FILE, true,
               "SDFiles with molfile for tagging as well " +
               InputFieldName.TAG_NAME.toString() + " and " +
               InputFieldName.SET_NAME.toString() + " fields." +
               " NOTE: SMARTS are read in first, then the SD file.");
      opt.setRequired( false );
      options.addOption( opt );

      opt = new Option( OPT_TAG_SETS, true,
               TAG_SETS_ALL +
               " or name1|name2|...; The given names must match the names" +
               " in the files disregarding the letter case" );
      opt.setRequired( true );
      options.addOption( opt );

      opt = new Option( OPT_TAG_INFO, true,
               OutputTag.getTagsAsText( "|" ) +
               " or if not specified all the listed ones.");
      opt.setRequired( false );
      options.addOption( opt );

      opt = new Option( OPT_OUTPUT_COUNT, false,
               "Add for each structure pattern a field containing the" +
               " occurrence of the given pattern in the molecule" +
               " (NOTE: if multiple patterns have the same tag name the" +
               " counts will be summed up)." );
      opt.setRequired( false );
      options.addOption( opt );

      opt = new Option( OPT_OUTPUT_EXISTS, false,
               "Add for each structure pattern a field containing 0 " +
               " if the pattern does not match and 1 if the pattern matches " +
               " (NOTE: if multiple patterns have the same tag name the" +
               " any match will yield a 1)." );
      opt.setRequired( false );
      options.addOption( opt );

      opt = new Option( OPT_ONLY_MATCHES, false,
               "Remove records from output unless they match at least one query pattern." );
      opt.setRequired( false );
      options.addOption( opt );
      
      opt = new Option( OPT_NO_MATCHES, false,
               "Remove records from output if they match any query pattern." );
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
      String molPatternFile = cmd.getOptionValue( OPT_QSDF_FILE );
      String tagSets = cmd.getOptionValue( OPT_TAG_SETS );
      String tagInfo = cmd.getOptionValue( OPT_TAG_INFO );
      boolean outputCount = cmd.hasOption( OPT_OUTPUT_COUNT );
      boolean outputExists = cmd.hasOption( OPT_OUTPUT_EXISTS );
      boolean onlyMatches = cmd.hasOption( OPT_ONLY_MATCHES );
      boolean noMatches = cmd.hasOption( OPT_NO_MATCHES );
      SDFStructureTagger tagger = new SDFStructureTagger( outFile );

      try
      {  tagger.prepare( smartsFile, molPatternFile, tagSets, tagInfo, outputCount, outputExists, onlyMatches, noMatches );
         tagger.run( inFile );
      } catch( IncorrectInputException iie )
      {  System.err.println( iie.toString() );
         exitWithHelp( options );
      } finally
      {  tagger.close();
      }
   }
}
