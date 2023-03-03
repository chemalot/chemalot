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

import java.beans.Statement;
import java.io.File;
import java.io.IOException;
import java.util.regex.Pattern;

import openeye.oechem.*;

import org.apache.commons.cli.*;
import org.jdom.JDOMException;

import com.aestel.Settings;
import com.aestel.io.dataAccess.*;
import com.aestel.io.dataAccess.Record;

/**
 * SDFExporter Exports the output of a sql statement to an sd file
 *
 * @author albertgo 2008
 *
 */
public class SDFSDFExporter
{  private static final Pattern NEWLinePattern = Pattern.compile("[\r\n]+");

   private final SQLStatement rowStmt;
   private final Selecter rowSel;
   private final String newLineReplacement;

   private int inCount = 0;
   private int nStruct = 0;
   private final long start;

   private final oemolithread in;
   private final oemolothread out;

   private final String[] tagNames;

   private final boolean printIfNoRecord;

   public static SDFSDFExporter createFromFile(
            String sqlFile, String rowStmtName,
            String inFile, String outFile, String [] tagNames,
            boolean printIfNoRecord, String newLineReplacement) throws IOException
   {  if( sqlFile == null ) sqlFile = "sql.xml";

      File sqlFILE = new File(sqlFile);
      if( ! sqlFILE.exists() )
      {  sqlFILE = new File(Settings.AESTEL_INSTALL_PATH + File.separatorChar
                           + "config" + File.separatorChar + "sdfExport",
                            sqlFile);
         if( ! sqlFILE.exists() )
            sqlFILE = new File(Settings.AESTEL_INSTALL_PATH + File.separatorChar
                     + "config" + File.separatorChar + "sdfExport",
                  sqlFile + ".xml" );


         if( ! sqlFILE.exists() )
            throw new IOException(sqlFile + " not found!");
      }
      SQLStatement stmt = SQLStatement.createFromFile(sqlFILE, rowStmtName);

      return new SDFSDFExporter(stmt, inFile, outFile, tagNames, printIfNoRecord, newLineReplacement);
   }


   public static SDFSDFExporter createFromStatementStr(
            String selectStr, String paramTypes,
            String inFile, String outFile, String [] tagNames,
            boolean printIfNoRecord, String newLineReplacement)
   {  String[] typeStr = paramTypes.split("\\s+");
      DADataType[] pTypes = new DADataType[typeStr.length];
      int i = 0;
      for(String pStr : typeStr)
         pTypes[i++] = Enum.valueOf(DADataType.class, pStr.toUpperCase());

      i = 0;
      DAParameterUse[] usg = new DAParameterUse[pTypes.length];
      for(String pStr : typeStr)
         usg[i++] = DAParameterUse.IN;

      SQLStatement stmt = SQLStatement.createStatement(selectStr, pTypes, usg );

      return new SDFSDFExporter(stmt, inFile, outFile, tagNames, printIfNoRecord, newLineReplacement);
   }


   private SDFSDFExporter(SQLStatement stmt,
            String inFile, String outFile, String [] tagNames,
            boolean printIfNoRecord, String newLineReplacement)
   {  start = System.currentTimeMillis();
      this.newLineReplacement = newLineReplacement;
      this.rowStmt = stmt;
      rowSel = Selecter.factory( rowStmt );

      in = new oemolithread(inFile);
      out = new oemolothread(outFile);

      this.tagNames = tagNames;
      this.printIfNoRecord = printIfNoRecord;
   }

   /*
    * SDFExporter -sqlFile sql.xmlFile -sqlName nameOfSQLInFile
    * outputFile [args]
    *
    * Will execute the sql command named nameOfSQLInFile in the xml file
    * sql.xml File. The sql is expected to return the rows with data to be
    * exported.
    *
    * sql fields are added to the input sdf file.
    */
   private static final String CMDLineText = "SDFSDFExport";
   private static final String HEADERText =
           "Add fields from sql query to SDF file. \n"
         + "Parameters for the sql statment are from tags in the input sdf file.\n"
         + "Use 'title' as tagName ofr column name to affect the molfiel title.\n\n";

   public static void main(String[] args) throws ParseException, JDOMException, IOException
   {  // create command line Options object
      Options options = new Options();
      Option opt = new Option("sqlFile", true, "sql-xml file");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("sqlName", true, "name of SQL element in xml file, Default sql.xml in 'sdfExport' config");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("selectStatement", true, "select statement to execute");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("paramTypes", true, "'|' separated list of parameter types to pass tostatment int,float,string,date");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("o", "out", true, "output file");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("i", "in", true, "input file, oe or .tab each record executes the query once. Use '.tab' to read from stdin");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("queryTags", true, "'|' separetaed list of tags whose values is passt to the sql.");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("newLineReplacement",true, "If given newlines in fields will be replaced by this string.");
      options.addOption(opt);

      opt = new Option("filterIfNoRecords", false, "If no rows are returned by the query that record is filtered out.");
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
      String selStr  = cmd.getOptionValue("selectStatement");
      String pTypes  = cmd.getOptionValue("paramTypes");
      String newLineReplacement = cmd.getOptionValue("newLineReplacement");
      boolean printIfNoRecord = ! cmd.hasOption("filterIfNoRecords");
      String[] tagStr = cmd.getOptionValue("queryTags").trim().split("\\|");

      try
      {
         SDFSDFExporter exporter = null;
         if(  (sqlFile != null && sqlFile.length() > 0) || (sqlName != null && sqlName.length() > 0))
         {  if( (selStr != null && selStr.length() > 0) || (pTypes != null && pTypes.length() > 0))
            {  System.err.println("sqlFile and sqlName may not be used with selectStatement and paramTypes");
               exitWithHelp(options);
            }

            exporter = createFromFile(sqlFile, sqlName, inFile, outFile, tagStr,
                                       printIfNoRecord, newLineReplacement);

         }else if( selStr == null || selStr.length() == 0 || pTypes == null || pTypes.length() == 0)
         {  System.err.println("sqlFile and sqlName or selectStatement and paramTypes must be given");
            exitWithHelp(options);

         }else
         {  exporter = createFromStatementStr(selStr, pTypes, inFile, outFile, tagStr,
                                       printIfNoRecord, newLineReplacement);
         }

         exporter.export();
         exporter.close();
      } catch (Exception e)
      {  e.printStackTrace();
         System.err.println();
         exitWithHelp(options);
      }

   }


   private void export() throws IOException
   {  OEGraphMol mol = new OEGraphMol();

      while(oechem.OEReadMolecule(in, mol))
      {  String[] param = getQueryValues(mol);
         export(mol, param);
      }

      mol.delete();
   }

   private boolean export(OEGraphMol mol, String[] args) throws IOException
   {  inCount++;

      if (!rowSel.select((Object[]) args))
      {  if( printIfNoRecord )
            oechem.OEWriteMolecule(out, mol);
         nStruct++;
         return false;
      }

      String[] fieldNames = rowSel.getFieldNames();

      while (rowSel.hasNext())
      {  Record sqlRec = rowSel.next();
         OEGraphMol mcopy = new OEGraphMol( mol );

         for (int i = 0; i < fieldNames.length; i++)
         {  String fld = sqlRec.getStrg(i).trim();
            if( fld.length() == 0) continue;

            if( newLineReplacement != null )
               fld = NEWLinePattern.matcher(fld).replaceAll(newLineReplacement);
            if( "title".equalsIgnoreCase(fieldNames[i]))
               mol.SetTitle(fld);
            else
               oechem.OESetSDData(mcopy, fieldNames[i], fld);
         }

         oechem.OEWriteMolecule(out, mcopy);
         mcopy.delete();

         nStruct++;
      }

      return true;
   }

   private String[] getQueryValues(OEGraphMol mol)
   {  String[] ret = new String[tagNames.length];
      for(int i=0; i<tagNames.length; i++)
      {  if( "title".equalsIgnoreCase(tagNames[i]) )
            ret[i] = mol.GetTitle();
         else
            ret[i] = oechem.OEGetSDData(mol, tagNames[i]);
      }
      return ret;
   }

   public void close()
   {  System.err.printf("SDFExporter: Input %d, exported %d structures in %dsec\n",
             inCount, nStruct,
            (System.currentTimeMillis() - start) / 1000);

      if(rowSel != null) rowSel.close();
      in.close();
      in.delete();
      out.close();
      out.delete();
   }


   private static void exitWithHelp(Options options)
   {  HelpFormatter formatter = new HelpFormatter();
      formatter.printHelp(80, CMDLineText, HEADERText, options, "");
      System.exit(1);
   }
}
