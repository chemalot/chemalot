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
package com.genentech.chemistry.tool;

import java.io.IOException;
import java.util.LinkedHashSet;
import java.util.Set;

import openeye.oechem.*;

import org.apache.commons.cli.*;
import org.jdom.JDOMException;

import com.aestel.Settings;
import com.aestel.chemistry.depict.DepictHelper;
import com.aestel.chemistry.depict.ImageType;
import com.genentech.oechem.tools.OETools;

public class SDF2HtmlTab
{
   /*
    * SDF2HtmlTab -in fName.sdf >t.xls
    *
    * Will read the sd file and create a tab separated file with columns for
    * each tag in the sd file. The sd file is read twice, first to find the tags.
    *  Then to create the html.
    */
   private static final String INTROText = "SDF2HtmlTab -in fName.sdf\n"
      +"Create html fiel with image tags using the genetech webserver for structures.\n"
      +"\n";

   private static final String BASEUrl = "http://research/";
   private static int IMGWidth  = 120;
   private static int IMGHeigth = 120;

   public static void main(String [] args) throws ParseException, JDOMException, IOException {
      long start = System.currentTimeMillis();

      // create command line Options object
      Options options = new Options();
      Option opt = new Option("in",true, "input sd file");
      opt.setRequired(true);
      options.addOption(opt);

      CommandLineParser parser = new BasicParser();
      CommandLine cmd = null;
      try {
         cmd = parser.parse( options, args);
      } catch(Exception e) {
         System.err.println(e.getMessage());
         exitWithHelp(options);
      }
      args = cmd.getArgs();

      String inFile = cmd.getOptionValue("in");

      args = cmd.getArgs();
      if(args.length > 0) {
         exitWithHelp(options);
      }

      oemolistream ifs;
      Set<String> tagSet = getTagSet(inFile);

      System.out.println("<html xmlns:v='urn:schemas-microsoft-com:vml' xmlns:o='urn:schemas-microsoft-com:office:office'>");
      System.out.println("<head>");
      System.out.println("<BASE href='"+BASEUrl+"'/>");
      System.out.println("<link href='/" + Settings.SERVLET_CONTEXT +
                         "/css/Aestel.css' rel='stylesheet' type='text/css'/>");

      System.out.println("<style type='text/css'>" );
      System.out.println("td.stru { width: "+(IMGWidth+2)+"px; height: "+(IMGHeigth+4)+"px; vertical-align: top; }");
      System.out.println("table.grid tr.first { border-top: 3px solid black; }");
      // for tables in tables
      System.out.println("table.grid table td { border: 0px; text-align: right;}");
      System.out.println("table.grid table td:first-child { text-align: left;}");
      System.out.println(
         "th.head { border-left: 1px solid black; border-bottom: 2px solid black;\n"
       + "          empty-cells: show; background-color: #6297ff; color: #000000;\n"
       + "          padding: 0em .3em 0em .3em; vertical-align: middle; }");
      System.out.println("</style>" );

      System.out.println("</head>" );
      System.out.println("<body>" );

      System.out.println("<table class='grid'><tr>" );
      System.out.println("<th class='head'>Structure</th>" );

      for(String tag : tagSet)
         System.out.println("<th class='head'>" + tag + "</th>");
      System.out.println("</tr>");

      OEGraphMol mol = new OEGraphMol();
      ifs = new oemolistream(inFile);
      int iCounter = 0;
      while(oechem.OEReadMolecule(ifs , mol) ) {
         iCounter++;
         System.out.println("<tr>");
         String smi = OETools.molToCanSmi(mol, true);
         String img = DepictHelper.DEFAULT.getExcelSmilesImageElement(
                                    BASEUrl, 120, 120, ImageType.PNG, smi, null);

         System.out.print(" <td class='stru'>");
         System.out.print(img);
         System.out.println("</td>");

         for(String tag : tagSet) {
            String val = oechem.OEGetSDData(mol, tag);

            System.out.print(" <td>");
            System.out.print(val);
            System.out.println("</td>");
         }

         System.out.println("</tr>");
      }

      System.out.println("</table></body></html>");

      System.err.printf("SDF2HtmlTab: Exported %d structures in %dsec\n",
               iCounter, (System.currentTimeMillis()-start)/1000);
   }

   private static Set<String> getTagSet(String inFile)
   {  oemolistream ifs = new oemolistream(inFile);

      OEGraphMol mol = new OEGraphMol();
      Set<String> tagSet = new LinkedHashSet<String>();

      while(oechem.OEReadMolecule(ifs , mol ) ) {
         OESDDataIter sdPairs = oechem.OEGetSDDataPairs(mol);
         while(sdPairs.hasNext()) {
            OESDDataPair sdPair = sdPairs.next();
            tagSet.add(sdPair.GetTag());
         }
      }
      ifs.close();
      ifs.delete();
      mol.delete();

      return tagSet;
   }

   private static void exitWithHelp(Options options) {
      HelpFormatter formatter = new HelpFormatter();
      formatter.printHelp( INTROText, options );
      System.exit(1);
   }
}

