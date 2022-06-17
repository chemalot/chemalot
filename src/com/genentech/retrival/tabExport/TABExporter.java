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
package com.genentech.retrival.tabExport;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.regex.Pattern;

import org.apache.commons.cli.*;
import org.jdom.JDOMException;

import com.aestel.io.dataAccess.DAError;
import com.aestel.io.dataAccess.Record;
import com.aestel.io.dataAccess.SQLStatement;
import com.aestel.io.dataAccess.Selecter;

/**
 * SDFExporter Exports the output of a sql statement to an sd file
 * @author albertgo 2008
 *
 */
/*
 * TABExporter -sqlFile sql.xmlFile -sqlName nameOfSQLInFile
 *             -o outputFile [args]
 *
 * Will execute the sql command named nameOfSQLInFile in the xml file sql.xmlFile.
 * The sql is expected to return the rows with data to be exported.
 *
 * The sql statement is executed and for each rows a record of the tab file is
 * containing the molfile and one tag for each of the other columns in the sql
 * statement. The column name is used as TAB header-name and the column values as
 * the row values.
 */
public class TABExporter
{
   private static final Pattern NEWLinePattern = Pattern.compile("[\r\n]+");

   private static final String INTROText = "TABExport\n"
      +"Export TAB file from given sql statement.\n"
      +"If additional parameter are given these are passed to the sql statement.\n";


   public static void main(String [] args) throws ParseException, JDOMException, IOException
   {  long start = System.currentTimeMillis();
      int nRows = 0;

      // create command line Options object
      Options options = new Options();
      Option opt = new Option("sqlFile",true, "sql-xml file");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("sqlName",true, "name of SQL element in xml file");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("o", "out", true, "output file");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("i", "in", true, "tab separated file with input arguments for query");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("newLineReplacement",true, "If given newlines in fields will be replaced by this string.");
      options.addOption(opt);

      opt = new Option("noHeader",false, "Do not output header line");
      options.addOption(opt);

      CommandLineParser parser = new BasicParser();
      CommandLine cmd = null;
      try
      {  cmd = parser.parse( options, args);
      } catch(Exception e)
      {  System.err.println(e.getMessage());
         exitWithHelp(options);
      }
      
      String outFile = cmd.getOptionValue("o");
      String inFile  = cmd.getOptionValue("i");
      String sqlFile = cmd.getOptionValue("sqlFile");
      String sqlName = cmd.getOptionValue("sqlName");
      String newLineReplacement = cmd.getOptionValue("newLineReplacement");

      args = cmd.getArgs();
      if( inFile != null && args.length > 0)
      {  System.err.println("-i does not work with additional arguments");
         exitWithHelp(options);
      }


      try {
         PrintStream out = System.out;
         if(outFile != null) out = new PrintStream(outFile);

         boolean printHeader = ! cmd.hasOption("noHeader");

         
         SQLStatement stmt = SQLStatement.createFromFile(new File(sqlFile), sqlName);
         
         if( inFile == null )
         {  nRows = exportByArgs(stmt, args, out, printHeader, newLineReplacement);
         }else
         {  BufferedReader in;
            if( ".tab".equalsIgnoreCase(inFile) )
               in = new BufferedReader(new InputStreamReader(System.in));
            else
               in = new BufferedReader(new FileReader(inFile));

            String line;
            while((line=in.readLine()) != null)
            {  args = line.split("\t");
            
               if( nRows > 0 ) printHeader = false;
               
               nRows += exportByArgs(stmt, args, out, printHeader, newLineReplacement);
            }
            in.close();
         }

      } catch (Exception e) {
         throw new Error(e);
      }finally {
         System.err.printf("TABExporter: Exported %d records in %dsec\n",
               nRows, (System.currentTimeMillis()-start)/1000);
      }
   }

   private static int exportByArgs(SQLStatement stmt, String[] args,
         PrintStream out, boolean printHeader, String newLineReplacement) throws DAError
   {  int nRows = 0;
   
      Object[] sqlArgs = args;
      if(stmt.getParamTypes().length != args.length)
      {  System.err.printf("\nWarining sql statement needs %d parameters but got only %d. Filling up with NULLs.\n",
            stmt.getParamTypes().length, args.length);
         sqlArgs = new Object[stmt.getParamTypes().length];
         System.arraycopy(args, 0, sqlArgs, 0, args.length);
      }

      Selecter sel = Selecter.factory(stmt);
      if(! sel.select(sqlArgs) )
      {  System.err.println("No rows returned!");
         System.exit(0);
      }

      String[] fieldNames = sel.getFieldNames();
      if( fieldNames.length == 0)
      {  System.err.println("Query did not return any columns");
         sel.close();
         return 0;
      }

      if( printHeader )
      {  StringBuilder sb = new StringBuilder(200);
         for(String f : fieldNames)
            sb.append(f).append('\t');
         if(sb.length() > 1) sb.setLength(sb.length()-1);  // chop last \t
         String header = sb.toString();

         out.println(header);
      }

      StringBuilder sb = new StringBuilder(200);
      while ( sel.hasNext() )
      {  Record sqlRec = sel.next();
         sb.setLength(0);

         for(int i=0; i<fieldNames.length; i++)
         {  String fld = sqlRec.getStrg(i);
            if( newLineReplacement != null )
               fld = NEWLinePattern.matcher(fld).replaceAll(newLineReplacement);

            sb.append(fld).append('\t');
         }

         if(sb.length() > 1) sb.setLength(sb.length()-1);  // chop last \t
         String row = sb.toString();

         out.println(row);

         nRows++;
      }
      sel.close();
      
      return nRows;
   }

   private static void exitWithHelp(Options options) {
      HelpFormatter formatter = new HelpFormatter();
      formatter.printHelp( INTROText, options );
      System.exit(1);
   }
}
