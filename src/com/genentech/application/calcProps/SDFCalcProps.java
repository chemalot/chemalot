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
package com.genentech.application.calcProps;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.jdom.Element;

import com.aestel.Settings;
import com.aestel.io.IOUtil;
import com.aestel.io.XMLUtil;
import com.aestel.utility.exception.IncorrectInputException;

/**
 * Calculates properties of small molecules using various command line programs.
 * Properties and their respective command line programs are specified in XML files.
 * Properties have dependencies and programs must be sorted in order of dependency.
 * For example, Solubility_Index depends on cLogD7.4 and Aromatic_Ring. These two 
 * properties must be calculated before calculating solubility index.
 *
 * @author JW Feng & Man-Ling Lee / Last updated on Oct 31, 2015
 * Copyright 2012-2015 Genentech
 */
public class SDFCalcProps
{
   /** 
    * Get a list of calculators based on the property requested by the user
    * 
    * getCalculators() calls itself in case the requiredCalculator list is not
    * empty.
    * 
    * @param prop
    *           Name of the requested properties
    * @param availCALCS
    *           is the full list of calculators defined by the XML files
    * @return  a set of calculators 
    */
   private static Set<Calculator> getCalculators(String prop, Set<Calculator> availCALCS)
   {
      Set <Calculator> myCalcs = new LinkedHashSet<Calculator>();
      for (Calculator calc : availCALCS){
         if (calc.getName().equals(prop)){
            myCalcs.add(calc);
            for (String dep : calc.getRequiredCalculators()){
               //recursively get all dependent calculators
               myCalcs.addAll(getCalculators(dep, availCALCS));
            }
         }
      }
      return myCalcs;
   }

   /**Generated a string of piped commands based on a set of calculators
    * Need to figure out dependencies,
    */
   private static String assembleCommands(Set<Calculator> calculators, boolean verbose, 
            boolean debug, String counterTag, Set<String> outputTags) {
      if (calculators.isEmpty()) {
         return null;
      }

      //generate a string of SD tags that are separated by | character
      StringBuilder tagsToKeep = assembleTags(outputTags);

      if (debug) { //extra info debugging
         System.err.println("++++++++++Calculators for sorting ++++++++++++++++");
         for (Calculator c: calculators) {
            Set<String> deps = c.getRequiredCalculators();
            System.err.print(c.getName() + ", " + c.getProgName() + " " + 
                     c.getProgOps() + ", number of required calculators: " + 
                     deps.size() + "(");
            for (String d : deps) 
               System.err.print(d + " ");
            System.err.println(")");
         }
      }

      //sort calculators by dependencies
      //cHLM required cLogD7.4, cLogD7.4 should be calculated prior to cHLM
      List <Calculator> sortedCalculators = sortByDependencies(
               new ArrayList<Calculator> (calculators), 0);

      if (debug) { //extra info for debugging
         System.err.println("++++++++++sorted Calculators++++++++++++++++");
         for (Calculator c: sortedCalculators) {
            Set<String> reqCalculators = c.getRequiredCalculators();
            System.err.print(c.getName() + ", " + c.getProgName() + " " + 
                     c.getProgOps() + ", number of required calculators: " + 
                     reqCalculators.size() + "(");
            for (String d : reqCalculators) 
               System.err.print(d + " ");
            System.err.println(")");
         }
      }

      if (debug) {//extra info for debugging
         System.err.println("++++++++++aggregated Calculators++++++++++++++++");
         for (Calculator c: sortedCalculators) {
            Set<String> reqCalculators = c.getRequiredCalculators();
            System.err.print(c.getProgName() + " " + 
                     c.getProgOps() + ", number of required calculators: " + 
                     reqCalculators.size() + "(");
            for (String d : reqCalculators) 
               System.err.print(d + " ");
            System.err.println(")");
         }
      }

      //chain together the calculators to produce a single command line
      String calcCommands =  getCommands(sortedCalculators) +
                             " | sdf2Tab.csh -in .sdf -tags \"" + counterTag + "|" + tagsToKeep + "\"" ;
      return calcCommands;
   }


   /**
    * Get the list of SD tags that a property produces
    * also gets the tags that are produced by calculators this property depends on aka "required tags"
    * "required tags" are defined in the keepRequiredCalculators field on the XML files
   */
   private static TreeSet<String> getOutputFields(String prop, 
            Set<Calculator> calculators, boolean verbose)
   {
      TreeSet <String> tags = new TreeSet<String>(String.CASE_INSENSITIVE_ORDER);
      for (Calculator calc : calculators){
         if(calc.getName().equals(prop)){
            //add tags for calculator
            tags.addAll( calc.getOutputFields() );
            if (verbose) {
                tags.addAll( calc.getVerboseFields() );
            }
            //get tags for keepRequiredCalculators, recursively
            Set<String> requiredCalcs = calc.getKeepRequiredCalculators();
            for(String c : requiredCalcs ){
               tags.addAll(getOutputFields(c, calculators, verbose));
            }
         }
      }
      return tags;
   }

   /**
    * Generate a string of piped commands
    */
   private static String getCommands(List<Calculator> sortedCalcs)
   {
      String command=null;
      for(Calculator calc : sortedCalcs){
         if (calc.getProgName().length() > 0){
            if (command ==null ){//first time
               command = calc.getProgName() + " " + calc.getProgOps();
            }else{
               command = command + " | " + calc.getProgName() + " " + calc.getProgOps();
            }
         }
      }
      if (command ==null ){
         command = "tee";
      }
      return command;
   }
      

   /**
    * Sort commands by dependencies ex. Solubility_Index requires cLogD7.4
    * Need to calculate cLogD7.4 before calculating solubility_index
    * cLogD7.4 command line need to appear before solubility_index command line
    */
   private static List<Calculator> sortByDependencies(
            ArrayList<Calculator> calculators, int noCalculatorSizeChangeCount)
   {
      int oldCalculatorsSize = calculators.size();
      List<Calculator> sorted = new ArrayList<Calculator>();
      
      if (! calculators.isEmpty()){
         Calculator calc = calculators.remove(0); //get first element in list

         Set<String> reqCalculators = calc.getRequiredCalculators();
         if (reqCalculators.size() == 0){
            // no dependencies, add to beginning
            sorted.add(0, calc);
         }else { //there are dependencies
            // are any dependencies left in the list of calculators to be sorted
            if (anyDependenciesInList(reqCalculators, calculators)){
               calculators.add(calc); //add calc back to the end of the list to be sorted later
            }else {
               //they must be in the sorted list, add calc to the end of sorted list
               sorted.add(calc); //append to end of sorted calculators
            }
         }
      }
      if( calculators.size() == oldCalculatorsSize )
         noCalculatorSizeChangeCount = noCalculatorSizeChangeCount + 1;
      else
         noCalculatorSizeChangeCount = 0;
      
      /*If the number of calculators in the list has not going down within
        calculators.size() times*/
      if( noCalculatorSizeChangeCount == calculators.size() 
       && calculators.size() > 0 )
      {  StringBuffer calculatorText = new StringBuffer();
         for( Calculator calc : calculators )
            calculatorText = calculatorText.append( calc.getName() ).append( " " );
         throw new Error( "There is a circular dependencies amongst following calculators: "
                  + calculatorText.substring( 0, calculatorText.length() ) );
      }

      //recursively sort remaining calculators
      if (calculators.size() > 0) {
         //append rest to sorted
         sorted.addAll(sortByDependencies(calculators, noCalculatorSizeChangeCount)); 
      }
      return sorted;
   }


   /**
    * Check if any of the required calculators in the list of calculators
    * 
    * @param requiredCalculators
    *          Set of required calculators retrieved from a calculator
    * @param calculators
    *          Assembled calculators for compiling the command text
    * @return  true if one required calculator exists in the assemble calculator
    *          list
    */
   private static boolean anyDependenciesInList( 
              Set<String> requiredCalculators, ArrayList<Calculator> calculators )
   {  boolean found = false;
      for (String d : requiredCalculators){
         for (Calculator calc : calculators){
            if (d.equals(calc.getName()))
               return true;
         }
      }
      return found;
   }


   private static String calculate(String[] props, boolean predictTautomer,
            boolean dontFilter, boolean verbose, boolean debug,
            boolean printOnly, boolean addMolIndex,
            Set<Calculator> availCALCS, String inFile, String outFile)
                     throws IOException, InterruptedException
   {
      String counterTag = "___sdfCalcProps_counter___";
      String savedTitleTag = "___sdfCalcProps_saved_title___";
      String tempFileRoot =  "$TMPDIR/sdfCalcProps.$$." + System.currentTimeMillis() ;
      String tempOrigFileName = tempFileRoot + ".orig.sdf";
      String filteredFileName = tempFileRoot + ".filtered.sdf";

      Set<Calculator> calculators = new LinkedHashSet<Calculator>();
      
      //Properties that depend on ionization state of the molecule ex. charge
      Set<Calculator> ionizedCalculators = new LinkedHashSet<Calculator>();
      //Properties that depend on the neutral molecule. ex. MW
      Set<Calculator> neutralCalculators = new LinkedHashSet<Calculator>();

      for (String prop : props){
         //getCalculators takes care of getting all dependent calculators, recursively
         Set<Calculator> myCalcs = getCalculators(prop, availCALCS);

         //adding to a set, does not contain duplicates calculators
         calculators.addAll(myCalcs);
      }

      //divide calculators into those that require ionization and those that don't
      for (Calculator calc : calculators) {
         if (calc.requiresIonization()) {
            ionizedCalculators.add(calc);
         }else {
            neutralCalculators.add(calc);
         }
      }

      //get the set of SD tags that will be produced, each property produces a set of SD tags
      //using TreeSet to keep the SD tags in alphabetical order
      TreeSet <String> ionizedOutputTags = new TreeSet<String>();
      for (Calculator p: ionizedCalculators){
          ionizedOutputTags.addAll(getOutputFields(p.getName(), ionizedCalculators, verbose));
      }

      TreeSet <String> neutralOutputTags = new TreeSet<String>();
      for (Calculator p: neutralCalculators){
          neutralOutputTags.addAll(getOutputFields(p.getName(), neutralCalculators, verbose));
      }

      // The list of SD tags that will be produced
      // this set should be a union of tags from ionized and neutral tags
      // I guess I could just merge the treesets from ionized and neutral tags
      TreeSet <String> allOutputTags = new TreeSet<String>(String.CASE_INSENSITIVE_ORDER);
      for (String prop : props){
         allOutputTags.addAll(getOutputFields(prop, calculators, verbose));
       }


      if (debug)
      {
         System.err.println("=============================================");
         System.err.println("The following properties will be calculated.");
         printProperties(calculators, true);

         System.err.println("The following tags will be produced.");
         for(String t :allOutputTags){
            System.err.print(t + " ");
         }
         System.err.println();
      }

      //The following tags will be produced, exit after printing
      if (printOnly)
      {
         StringBuilder outputTags = assembleTags(allOutputTags);
         System.out.println("echo '" + outputTags +"'");
         System.exit(0);
      }

      // special properties that dictate how molecules are preprocessed
      Calculator tautomerCalculator = null;
      Calculator filterCalculator =null;
      Calculator ionizeCalculator =null;
      for (Calculator calc : availCALCS){
         if (calc.getName().equals("predictTautomer")){
            tautomerCalculator = calc;
         }
         if (calc.getName().equals("filter")){
            filterCalculator = calc;
         }
         if (calc.getName().equals("ionize")){
            ionizeCalculator = calc;
         }
      }

      
     // assemble the command line base on the properties that were requested,
     // this is the most complicated part of this program "assembleCommands"

     //get a string of piped commands for calculating properties that depend on ionization state
      ionizedCalculators = consolidateByAggregationId( ionizedCalculators );
      String ionizedCommand = assembleCommands(ionizedCalculators, verbose, debug, counterTag, ionizedOutputTags);
      if (ionizedCommand != null ) {
         //prepend command to generated ionized molecules
         ionizedCommand = ionizeCalculator.getProgName() + " " + ionizeCalculator.getProgOps() + " | "  + ionizedCommand;
      }

     //get a string of piped commands for calculating properties on the neutral molecule
      neutralCalculators = consolidateByAggregationId( neutralCalculators );
      String neutralCommand = assembleCommands(neutralCalculators, verbose, debug, counterTag, allOutputTags);

      //save a temp file that contains a unique identifier
      //run filter to get rid of "bad" molecules
      String command = "sdfTagTool.csh -copy TITLE=" + savedTitleTag + " -addCounter -counterTag " + counterTag +
               " -title " + counterTag + " -in " + inFile + " -out .sdf | tee " + tempOrigFileName + " | " +
               filterCalculator.getProgName() + " " + filterCalculator.getProgOps() +
               " | sdfTagTool.csh -in .sdf -out .sdf -keep " + counterTag;

      //different command if not filtering
      if (dontFilter) {
         command = "sdfTagTool.csh -copy TITLE=" + savedTitleTag + " -addCounter -counterTag " + counterTag +
                  " -title " + counterTag + " -in " + inFile + " -out .sdf | tee " + tempOrigFileName +
                  " | sdfTagTool.csh -in .sdf -out .sdf -keep " + counterTag;
      }

      String cleanUpCommand = "sdfTagTool.csh -in .sdf -title " + savedTitleTag + " -out .sdf " +
               " | sdfTagTool.csh -in .sdf -remove \"" + counterTag + "|" + savedTitleTag + "\" -out " + outFile;

      //command to create a Mol_Index tag for each molecule
      if (addMolIndex) {
             cleanUpCommand = "sdfTagTool.csh -in .sdf -title " + savedTitleTag + " -out .sdf " +
               " -format 'Mol_Index=Mol_{" + counterTag + "}'" +
               " | sdfTagTool.csh -in .sdf -remove \"" + counterTag + "|" + savedTitleTag + "\" -out " + outFile;
      }

      //predict tautomer
      if (predictTautomer) {
         command = command + " | " + tautomerCalculator.getProgName() + " " +
                  tautomerCalculator.getProgOps() +  " > " +  filteredFileName;
      }else {
         command = command  +  " > " +  filteredFileName;
      }

//      command = command + "; cat " + filteredFileName;
      command = command + "; echo NCCO | sdfTagTool.csh -in .smi -out .sdf >> " + filteredFileName + "; cat " + filteredFileName;

      if (ionizedCommand != null && neutralCommand !=null) {
      // merge temp file with the two tab files
         String mergeCommand =
               "sdfTabMerger.csh -outAll -addEmptyValues -sdf " + tempOrigFileName + " -tab -  -mergeTag " + counterTag + " -mergeCol " + counterTag + " -out .sdf | " +
                cleanUpCommand;
         command = command + " | " + ionizedCommand +
               " | sdfTabMerger.csh -outAll -addEmptyValues -tab - -sdf " + filteredFileName + " -mergeTag " + counterTag + " -mergeCol " + counterTag + "  -out .sdf | " +
                  neutralCommand +  " | " + mergeCommand;
      }else if (ionizedCommand == null && neutralCommand !=null) {
         String mergeCommand = "sdfTabMerger.csh -outAll -addEmptyValues -sdf " + tempOrigFileName + " -tab - -mergeTag " +
               counterTag + " -mergeCol " + counterTag + " -out .sdf | " + cleanUpCommand;

         command = command + " | " + neutralCommand + " | " + mergeCommand;
      }else if (ionizedCommand != null && neutralCommand ==null) {
         String mergeCommand = "sdfTabMerger.csh -outAll -addEmptyValues -sdf " + tempOrigFileName + " -tab - -mergeTag " +
               counterTag + " -mergeCol " + counterTag + " -out .sdf | " + cleanUpCommand;

         command = command + " | " + ionizedCommand + " | " + mergeCommand;
      }

      if (debug) {
         System.err.println("ionized command:\n" + ionizedCommand);
         System.err.println("neutral command:\n" + neutralCommand);
         System.err.println("Command to be executed:\n" + command);
      }

      return (command);
  }

   /**
    * Consolidate calculators with the same aggregation id into one calculator
    * with the given aggregation id as calculator name.
    * 
    * The new calculator contains the combined progOptions, requiredCalculators,
    * keepRequiredCalculators, outputFields, and verboseFields. The list of the 
    * original calculators are kept in the 
    * 
    * @param oldCalcSet
    *          The calculator set containing calculators with aggregation IDs
    * @return the consolidated calculator set in which the calculators with the 
    *          same aggregation ID are replaced by the corresponding aggregation 
    *          calculator. The requiredCalculators and keepRequiredCalculators
    *          lists in calculators containing the names of calculators with a
    *          aggregation ID are replaced with the name of their respective 
    *          aggregation calculators
    */
   private static Set<Calculator> consolidateByAggregationId( Set<Calculator> oldCalcSet )
   {  
      //Resulting calculator list after the aggregation
      Set<Calculator> newCalcSet = new LinkedHashSet<Calculator>();
      
      //List of calculators replacing the calculators with aggregation ID
      Map<String,Calculator> aggCalcMap = new HashMap<String,Calculator>();
      
      for( Calculator oldCalc : oldCalcSet )
      {  
         String oldAggId = oldCalc.getProgAggregateID();
         String oldCalcName = oldCalc.getName();
         
         if( oldAggId == null || oldAggId.length() == 0 )
         {  newCalcSet.add( oldCalc );
            continue;
         }
         
         Calculator aggCalc = aggCalcMap.get( oldAggId );
         if( aggCalc == null )
         {  
            oldCalc.setName( "AggCalc_" + oldAggId );
            oldCalc.addOriginalCalculator( oldCalcName );
            oldCalc.setHelpText( "Aggregation of " + oldCalcName );
            newCalcSet.add( oldCalc );
            aggCalcMap.put( oldAggId, oldCalc );
            
         } else
         {
            if( !aggCalc.getProgName().equalsIgnoreCase( oldCalc.getProgName() ) )
            {  throw new Error( "Expect program name \"" + aggCalc.getProgName() + 
                        "\" for calculators with aggregation Id \"" + oldAggId +
                        " but " + oldCalcName + " has the program name " + 
                        oldCalc.getProgName() + "\"." );
            }
            aggCalc.addOriginalCalculator( oldCalcName );
            aggCalc.setHelpText( aggCalc.getHelpText() + ", " + oldCalcName );
            aggCalc.setProgOps( aggCalc.getProgOps() + " " + oldCalc.getProgOps() );
            
            aggCalc.addRequiredCalculators( oldCalc.getRequiredCalculators() );
            aggCalc.addKeepRequiredCalculators( oldCalc.getKeepRequiredCalculators() );
            
            aggCalc.addOutputFields( oldCalc.getOutputFields() );
            aggCalc.addVerboseFields( oldCalc.getVerboseFields() );
         }
      }
      
      /*
       * For all calculators in newCalcSet, replace the calculator names in 
       * requiredCalculators and keepRequiredCalculatorsn lists with the names 
       * of the corresponding aggregation calculators
       */
      for( Calculator newCalc : newCalcSet )
      {
         Set<String> newReqCalcs     = newCalc.getRequiredCalculators();
         Set<String> newKeepReqCalcs = newCalc.getKeepRequiredCalculators();
         
         // Loop over the consolidation calculators
         Iterator<Calculator> aggCalcIterator = aggCalcMap.values().iterator();
         while( aggCalcIterator.hasNext() )
         {
            Calculator aggCalc = aggCalcIterator.next();
            String newAggCalcName = aggCalc.getName();
            
            /*
             * Remove the names of the aggregated (original) calculator from the 
             * requiredCalculators and keepRequiredCalculators lists. However,
             * only add the name of the current aggregation calculator newCalc
             * if removeAll() indicates change of lists and newCalc is not referring
             * to the same calculator instance as aggCalc. The latter condition
             * is essential for prevent circular dependencies. For exampled RO5
             * has the same aggregation ID as have some of its required calculators
             */
            if( newReqCalcs.removeAll( aggCalc.getOriginalCalculators() ) )
            {  if( ! newCalc.getName().equals( newAggCalcName ) )
                  newReqCalcs.add( newAggCalcName );
            }
            if( newKeepReqCalcs.removeAll( aggCalc.getOriginalCalculators() ) )
            {  if( ! newCalc.getName().equals( newAggCalcName ) )
                  newKeepReqCalcs.add( newAggCalcName );
            }
         }
      }
      return newCalcSet;
   }

   public static void main(String args[]){
      String usage = "sdfCalcProps [options] <list of space separated properties>\n";

      // create Options object
      Options options = new Options();
      // add  options
      options.addOption("h", false, "print help message.");
      options.addOption("help", false, "print help message.");
      options.addOption("in", true, "inFile in OE formats: Ex: a.sdf or .sdf");
      options.addOption("out", true, "outputfile in OE formats. Ex: a.sdf or .sdf ");
      options.addOption("useExp", false, "Use experimental values, if available, in property calculations. False by default.");
      
      Option setEnvVarOpt = new Option("setEnvVar", false, "Specify a name=value environment variable that can be used by wrapped programs.  e.g. MOKA_DIR=/mydir/moka/ or MOKA_MODEL= to unset.  Multiple env vars should be set with multiple -setEnvVar options.");
      setEnvVarOpt.setArgs(50);    // Max number of env vars to set.  Arbitrary, but needs to be > 1
      options.addOption(setEnvVarOpt);
      
      options.addOption("addMolIndex", false, "Creates a sd Tag called Mol_Index for each molecule where the values are 1 .. N");
      options.addOption("predictTautomer", false, "Run tauthor to predict the most likely tautomer and neutralize it. False by default.");
      options.addOption("print", false, "Do not run calculation, just print out the list of SD tags that will be produced. False by default");
      options.addOption("dontFilter", false, "Do not run filter to remove \"bad\" molecules. False by default.");
      options.addOption("propertiesFile", true, "Use this file instead of the default file properties.xml; " +
                        "it must be in the same folder as the default file, i.e. $AESTEL_DIR/config/properties");
      options.addOption("verbose", false, "Output verbose SD tags for each property. False by default.");
      options.addOption("debug", false, "Create a debug output SD file, not fully implemented. False by default.");
      options.addOption("showAll", false, "Print help for all properties, including non-public ones.");
      CommandLineParser parser = new PosixParser();

      try
      {
         CommandLine cmd = parser.parse(options, args);
         if (cmd.hasOption("h") || cmd.hasOption("help") ) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp( usage, options );
         }

         boolean useExp =false;
         if (cmd.hasOption("useExp") ) {
            useExp = true;
         }
         
         boolean showAll =false;
         if (cmd.hasOption("showAll") ) {
            showAll = true;
         }

         boolean predictTautomer=false;
         if (cmd.hasOption("predictTautomer")) {
            predictTautomer = true;
         }

         boolean dontFilter=false;
         if (cmd.hasOption("dontFilter")) {
            dontFilter = true;
         }

         boolean print=false;
         if (cmd.hasOption("print")) {
            print = true;
         }
         boolean verbose=false;
         if (cmd.hasOption("verbose")) {
            verbose = true;
         }

         boolean debug=false;
         if (cmd.hasOption("debug")) {
            debug = true;
         }

         boolean addMolIndex=false;
         if (cmd.hasOption("addMolIndex")) {
            addMolIndex = true;
         }

         String setEnvVarCommand = "";
         String[] envVarNameValuePairs = null;
         if (cmd.hasOption("setEnvVar") ) {
            envVarNameValuePairs = cmd.getOptionValues ("setEnvVar");
            for (String envVarNameValuePair : envVarNameValuePairs) {
               String[] envVarNameValue = envVarNameValuePair.split("=", 2);
               if (envVarNameValue == null || envVarNameValue.length != 2)
               {
                  System.err.println("setEnvVar syntax should be -setEnvVar myVar=myValue.  Multiple env vars should be set with multiple -setEnvVar options.");
                  exitWithHelp(usage, options);                  
               }              
               
               System.err.println("Setting custom environment variable: " + envVarNameValue[0] + " = " + envVarNameValue[1]);
               
               // Assuming cshell
               setEnvVarCommand += "setenv " + envVarNameValue[0] + " " + envVarNameValue[1] + ";";
            }
         }

         String filename = "properties.xml";
         String propertiesDir = Settings.AESTEL_INSTALL_PATH + "/config/properties";
         
         String propertiesFile = cmd.getOptionValue("propertiesFile");
         if (propertiesFile != null && propertiesFile.length() > 0) {
            if( ! propertiesFile.endsWith(".xml") )
            {
               System.err.println(propertiesFile + ": Properties files need to have xml file extension.");
               System.exit(1);
            }
            filename = propertiesFile;
         }
         File file = new File(propertiesDir + "/" + filename);
         if( ! file.exists() )
         {
            System.err.println(filename + ": Not found in " + propertiesDir);
            System.exit(1);
         }
               
         URL url = IOUtil.getConfigUrl(filename, propertiesDir,
                  "/com/genentech/application/calcProps", false);
         Element root = XMLUtil.getRootElement(url, true);
         
         //get set of available calculators as defined in properties.xml
         Set <Calculator> availCALCS = readXML(root, useExp);

         if (cmd.getArgList().isEmpty()){
            System.err.println("Enter at least one property");
            printProperties(availCALCS, showAll);
            exitWithHelp(usage, options);
         }

         if (showAll) {
            printProperties(availCALCS, showAll);
            exitWithHelp(usage, options);
         }

         String inFile = cmd.getOptionValue("in");
         if (inFile == null) {
            exitWithHelp(usage, options);
         }

         String outFile = cmd.getOptionValue("out");
         if (outFile == null) {
            exitWithHelp(usage, options);
         }

         //make sure requested args are valid
         if (containInvalidProps(cmd.getArgs(), availCALCS) ){
            printProperties(availCALCS, showAll);
            exitWithHelp(usage, options);
         }

         if (debug) {
            System.err.println("The following properties were requested:");
            for( String p : cmd.getArgs()){
               System.err.print(p+ " ");
            }
            System.err.println();
         }

         String[] props = cmd.getArgs();
         String command = calculate(props, predictTautomer, dontFilter, verbose, debug, print, addMolIndex,
                  availCALCS, inFile, outFile);
         
         // Prepend any environment variables 
         command = setEnvVarCommand + command;
         
         System.out.println(command);
                  

      } catch (ParseException e)
      {  // TODO print explanation
         throw new Error(e);
      } catch (IncorrectInputException e)
      {
         e.printStackTrace();
      } catch (IOException e)
      {
         throw new Error(e);
      } catch (InterruptedException e)
      {
         throw new Error(e);
      }
   }
   
   private static boolean containInvalidProps(String[] props, Set<Calculator> calculators)
   {
      boolean invalidProp = false;

      for (String prop : props){
         boolean found = false;
         for (Calculator calc : calculators){
            //if (calc.getName().equals(prop) && calc.isPublic()){
            if (calc.getName().equals(prop) ){
               found =true;
            }
         }
         if (found == false) {
            invalidProp = true;
            System.err.println(prop + " is not a valid property");
         }
      }

      return invalidProp;
   }

   public static Set<Calculator> readXML(Element root, boolean useExp){

      @SuppressWarnings("unchecked")
      List<Element> list = root.getChildren("property");
      Set<Calculator> calculators = new LinkedHashSet <Calculator>();

      for (int i = 0; i < list.size(); i++) {

         Element node = list.get(i);
         Calculator calc = new Calculator(node, useExp);
         calculators.add(calc);
      }

      //check to make sure all the required are valid
      for (Calculator calc : calculators){
         for (String reqCalc : calc.getRequiredCalculators()) {
            boolean found = false;
            for (Calculator calc2 : calculators){
               if (calc2.getName().equals(reqCalc)){
                  found=true;
               }
            }
            if (found == false){
               System.err.println("Invalid XML definition, " + reqCalc + " is a required calculator for " + calc.getName() +
                        " but " +  reqCalc + " is not a defined calculator.");
               System.exit(1);
            }
         }
      }

      //check to make sure all the required calculators are valid
      for (Calculator calc : calculators){
         for (String keepReqCalc : calc.getKeepRequiredCalculators()) {
            boolean found = false;
            for (Calculator calc2 : calculators){
               if (calc2.getName().equals(keepReqCalc)){
                  found=true;
               }
            }
            if (found == false){
               System.err.println("Invalid XML definition, " + keepReqCalc + " is a keep required calculator for " + calc.getName() +
                        " but " + keepReqCalc +" is not a defined calculator.");
               System.exit(1);
            }
         }
      }

      return calculators;
   }

   private static void printProperties(Set<Calculator> calculators, boolean showHidden)
   {
      //Print properties by alphabetical order
      TreeMap<String, String> sortedCalcs = new TreeMap<String, String>(String.CASE_INSENSITIVE_ORDER);
      for (Calculator calc : calculators){
         if (calc.isPublic()) {
            sortedCalcs.put(calc.getName(), calc.getHelpText() );
         }else if (showHidden == true) {//print non public props as well
            sortedCalcs.put(calc.getName(), calc.getHelpText() );
         }
      }

      for (String key: sortedCalcs.keySet()) {
         System.err.println(key + ":\t" + sortedCalcs.get(key));
      }
   }

   //generate a string of SD tags that are separated by | character
   private static StringBuilder assembleTags(Set <String> tags)
   {
      StringBuilder tagsToKeep = new StringBuilder();
      for (String tag : tags) {
         if (tagsToKeep.length() ==0) {
            tagsToKeep.append(tag);
         }else {
            tagsToKeep.append("|").append(tag);
         }
      }
      return tagsToKeep;
   }


   private static void exitWithHelp(String usage, Options options)
   {
      HelpFormatter formatter = new HelpFormatter();
      System.err.println(); System.err.println();

      formatter.printHelp( new PrintWriter(System.err, true), 100,  usage, "",  options, 2, 3, "");
      System.exit(1);
   }
}
