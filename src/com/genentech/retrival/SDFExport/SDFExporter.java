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
package com.genentech.retrival.SDFExport;

import java.io.*;
import java.util.regex.Pattern;

import openeye.oechem.OEGraphMol;

import org.apache.commons.cli.*;
import org.jdom.JDOMException;

import com.aestel.chemistry.depict.SmiTalk;
import com.aestel.io.dataAccess.*;
import com.aestel.utility.Message;
import com.aestel.utility.MessageList;
import com.aestel.utility.Message.Level;
import com.genentech.oechem.tools.OETools;
import com.genentech.struchk.oeStruchk.SaltRemover;

/**
 * SDFExporter Exports the output of a sql statement to an sd file
 *
 * @author albertgo 2008
 *
 */
public class SDFExporter
{  private static final Pattern NEWLinePattern = Pattern.compile("[\r\n]+");

   private final SQLStatement rowStmt;
   private final Selecter rowSel;
   private final Selecter molSel;
   private final PrintStream out;
   private final String molType;
   private final SmiTalk smiTalk;
   private final boolean ignoreMisMatches;
   private final boolean removeSalt;
   private final String newLineReplacement;
   private final MessageList msgs = new MessageList();

   private int inCount = 0;
   private int nStruct = 0;
   private SaltRemover saltRemover = null;
   private final long start;


   public static final String EMPTY_MOLFILE = "\n\n\n"
      + "  0  0  0  0  0  0  0  0  0  0999 V2000\n"
      + "M  END";

   private static final Pattern MOLFPattern = Pattern.compile("\\s\"?MOLFILE\"?(\\s|,)",
         Pattern.CASE_INSENSITIVE);
   private static final Pattern SMIPattern = Pattern.compile("\\s\"?SMILES\"?(\\s|,)",
         Pattern.CASE_INSENSITIVE);


   public SDFExporter(String sqlFile, String rowStmtName, String molStmtName,
                      String outFile, String newLineReplacement, boolean removeSalt,
                      boolean ignoreMismatches)
         throws IOException
   {  start = System.currentTimeMillis();
      this.newLineReplacement = newLineReplacement;

      File sqlFILE = new File(sqlFile);
      if( ! sqlFILE.exists() )
         throw new IOException(sqlFile + " not found!");

      rowStmt = SQLStatement.createFromFile(sqlFILE, rowStmtName);
      rowSel = Selecter.factory( rowStmt );
      rowSel.setFetchSize(200);
      this.ignoreMisMatches = ignoreMismatches;
      this.removeSalt = removeSalt;
      if (removeSalt)
         saltRemover = SaltRemover.createDefault();

      SQLStatement molStmt = rowStmt;
      if( molStmtName != null )
      {  molStmt = SQLStatement.createFromFile(sqlFILE, molStmtName);
         molSel = Selecter.factory( molStmt );

      }else
      {  molSel = null;
      }

      if( MOLFPattern.matcher(molStmt.toString()).find() )
      {  molType = "MOLFILE";
         smiTalk = null;

      }else if( SMIPattern.matcher(molStmt.toString()).find() )
      {  molType = "SMILES";
         smiTalk = new SmiTalk();
      }else
      {  throw new Error("Could not identify molecule format:\n" + molStmt.toString());
      }

      if (outFile == null || outFile.equalsIgnoreCase(".sdf"))
         out = System.out;
      else
         out = new PrintStream(outFile);
   }

   /*
    * SDFExporter -sqlFile sql.xmlFile -sqlName nameOfSQLInFile -removeSalt -o
    * outputFile [args]
    *
    * Will execute the sql command named nameOfSQLInFile in the xml file
    * sql.xmlFile. The sql is expected to return the rows with data to be
    * exported. The first column is expected to be called "MOLFILE" or "SMILES"
    * and contain the molfile/smiles to be used in the sd File. SmiTalk will be
    * used to convert a smiles.
    *
    * If removeSalt is given the molecule will be striped of any salts it
    * contains.
    *
    * The sql statement is executed and for each rows a record of the sd file is
    * containing the molfile and one tag for each of the other columns in the
    * sql statement. The column name is used as sdf- tag-name and the column
    * value as the sd field value. This is done using simple string
    * concatenation.
    */
   private static final String CMDLineText = "SDFExport";
   private static final String HEADERText =
           "Export SDF file with gne structures. NO OE-FORMATS SUPPORTED!\n"
         + "Parameters for the sql statment may be given in two ways:\n"
         +"   - as additinal arguments on the command line (one time execution)\n"
         +"   - as tab separated values in -inFile (one execution per in-line\n"
         + "The first column of the sql query must be a molfile or smiles and must"
         + "be called MOLFILE or SMILES or\n"
         + "you must specify the 'molSqlName' which is a select statement fetching a single"
         + "molfile/smiles column by the value of the first column in 'sqlName'"
         + "this column must be calles MOLFILE or SMILES.\n"
         + "Specifying a simple molfile query separatly with molSqlName can have "
         + "significant performance benefits for queries involving big sorts.";

   public static void main(String[] args) throws ParseException, JDOMException, IOException
   {  // create command line Options object
      Options options = new Options();
      Option opt = new Option("sqlFile", true, "sql-xml file");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("sqlName", true, "name of SQL element in xml file");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("molSqlName", true,
            "name of SQL element in xml file returning molfile for first column in sqlName");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("removeSalt", false,
            "remove any known salts from structure.");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("o", "out", true, "output file");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("i", "in", true, "input file, tab separated, each line executes the query once. Use '.tab' to read from stdin");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("newLineReplacement",true, "If given newlines in fields will be replaced by this string.");
      options.addOption(opt);

      opt = new Option("ignoreMismatches",false, "If not given each input must return at least one hit hits.");
      options.addOption(opt);


      CommandLineParser parser = new BasicParser();
      CommandLine cmd = null;
      try
      {  cmd = parser.parse(options, args);
      } catch (Exception e)
      {  System.err.println(e.getMessage());
         exitWithHelp(options);
      }

      String outFile = cmd.getOptionValue("o");
      String inFile  = cmd.getOptionValue("i");
      String sqlFile = cmd.getOptionValue("sqlFile");
      String sqlName = cmd.getOptionValue("sqlName");
      String molSqlName = cmd.getOptionValue("molSqlName");
      String newLineReplacement = cmd.getOptionValue("newLineReplacement");

      boolean removeSalt = cmd.hasOption("removeSalt");
      boolean ignoreMismatches = cmd.hasOption("ignoreMismatches");

      SDFExporter exporter = null;
      try
      {  exporter = new SDFExporter(sqlFile, sqlName, molSqlName, outFile,
                                    newLineReplacement, removeSalt, ignoreMismatches);
      } catch (Exception e)
      {  e.printStackTrace();
         System.err.println();
         exitWithHelp(options);
      }

      args = cmd.getArgs();

      if( inFile != null && args.length > 0 )
      {  System.err.println("inFile and arguments may not be mixed!\n");
         exitWithHelp(options);
      }

      if( inFile == null )
      {  if (exporter.getParamTypes().length != args.length)
         {  System.err.printf(
                  "sql statement (%s) needs %d parameters but got only %d.\n", sqlName,
                        exporter.getParamTypes().length, args.length);
            exitWithHelp(options);
         }

         exporter.export(args);
      }else
      {  BufferedReader in;
         if( ".tab".equalsIgnoreCase(inFile) )
            in = new BufferedReader(new InputStreamReader(System.in));
         else
            in = new BufferedReader(new FileReader(inFile));

         exporter.export(in);
         in.close();
      }

      exporter.close();
   }


   private void export(BufferedReader in) throws IOException
   {  String line;
      while( (line=in.readLine()) != null )
      {  String[] param = line.split("\t");
         if( !export(param) && ! ignoreMisMatches )
         {  break;
         }

         for(Message m : msgs.getMessages())
            System.err.println(m.getText());

         msgs.clear();
      }
   }

   private boolean export(String[] args) throws IOException
   {  inCount++;
      if (!rowSel.select((Object[]) args))
      {  msgs.addMessage(new Message(String.format("No rows returned for input %d!", inCount),
                         Level.WARNING, null));
         return false;
      }

      String[] fieldNames = rowSel.getFieldNames();
      if (fieldNames.length == 0)
      {  msgs.addMessage(new Message("Query did not return any columns.",
               Level.ERROR, null));
         return false;
      }

      OEGraphMol mol = new OEGraphMol();
      while (rowSel.hasNext())
      {  Record sqlRec = rowSel.next();
         String molStrOrID = sqlRec.getStrg(0);

         if (molStrOrID.length() == 0 || molStrOrID == null)
         {  molStrOrID = EMPTY_MOLFILE;

         } else
         {  if( molSel != null)
            {  if(!molSel.select(molStrOrID))
               {  msgs.addMessage(new Message("No Structure for: " + molStrOrID,
                        Level.WARNING, null));
                  molStrOrID = EMPTY_MOLFILE;
               }else
               {  Record molRec = molSel.next();
                  molStrOrID = molRec.getStrg(0);
               }
            }

            if ("SMILES".equals(molType))
               molStrOrID = smiTalk.smi2Mol(molStrOrID);
         }

         if (removeSalt)
            molStrOrID = removeSalts(molStrOrID);

         SDRecord sdRec = new SDRecord(molStrOrID, "");

         for (int i = 1; i < fieldNames.length; i++)
         {  String fld = sqlRec.getStrg(i);
            if( newLineReplacement != null )
               fld = NEWLinePattern.matcher(fld).replaceAll(newLineReplacement);
            sdRec.addData(fieldNames[i], fld);
         }

         out.printf("%s%s", sdRec, SDRecord.NEWLINE);

         nStruct++;
      }
      mol.delete();

      return true;
   }

   public void close()
   {  for(Message m : msgs)
         System.err.println(m.getText());
      msgs.clear();

      System.err.printf("SDFExporter: Exported %d structures in %dsec\n", nStruct,
            (System.currentTimeMillis() - start) / 1000);

      if(saltRemover != null ) saltRemover.delete();
      if(molSel != null) molSel.close();
      if(rowSel != null) rowSel.close();
      if(smiTalk != null) smiTalk.close();
   }

   private DADataType[] getParamTypes()
   {  return rowStmt.getParamTypes();
   }

   private String removeSalts(String molStr)
   {  OEGraphMol mol = new OEGraphMol();
      OETools.stringToMol(mol, molStr);
      MessageList mList = new MessageList();
      saltRemover.removeSalt(mol, msgs);

      molStr = OETools.molToString(mol);
      mol.delete();

      for (Message msg : mList.getMessages(Level.ERROR))
         System.err.println(msg);

      return molStr;
   }

   private static void exitWithHelp(Options options)
   {  HelpFormatter formatter = new HelpFormatter();
      formatter.printHelp(80, CMDLineText, HEADERText, options, "");
      System.exit(1);
   }
}
