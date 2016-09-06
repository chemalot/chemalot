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
package com.genentech.struchk;

import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.atomic.AtomicReference;
import java.util.regex.Pattern;

import javax.management.InvalidAttributeValueException;

import openeye.oechem.OEGraphMol;

import com.aestel.io.dataAccess.Record;
import com.aestel.io.dataAccess.SQLStatement;
import com.aestel.io.dataAccess.Selecter;
import com.genentech.oechem.tools.OETools;

/**
 * SDFExporter Exports the output of a sql statement to an sd file
 * @author albertgo 2008
 *
 */
public class SaltAdder {

//   private static final String INTROText = "SaltAdder\n"
//      +"\n";

   static final long RELAODTime = 3600 * 1000;

   private final OEGraphMol mol;
   private boolean isDeleted = false;

   public SaltAdder() {
      mol = new OEGraphMol();
   }
   
   /**
    * Add the structure of the salt (saltCode) to the molecule.
    * 
    * This will add multiple copies if the salt has a multiplicity.
    * 
    * @throws InvalidAttributeValueException if saltCode not defined
    */
   public void addSaltByCode(OEGraphMol mol, String saltCode) throws InvalidAttributeValueException {
      assert ! isDeleted;
      
      Salt salt = SaltMap.getSaltCodeMap().getSaltByCode(saltCode);
      if(salt == null)
         throw new InvalidAttributeValueException("Unknown Salt: " + saltCode);
      
      OETools.combineStructures(mol, salt.mol, salt.multiplicity);
   }
   
   /**
    * Add the structure of the salt (saltCode) to the molecule.
    * 
    * This will add multiple copies if the salt has a multiplicity.
    * 
    * @throws InvalidAttributeValueException if saltCode not defined
    */
   public String addSaltByCode(String molStr, String saltCode) throws InvalidAttributeValueException {
      assert ! isDeleted;
      
      mol.Clear();
      molStrToMol(mol, molStr);
      
      addSaltByCode(mol, saltCode);
      return OETools.molToString(mol);
   }

   /**
    * Add the structure of the salt (saltCode) to the molecule.
    * 
    * @param multiplicity number of copies to add to the molecule, the multiplicity 
    *    in the salt table is ignored.
    *    
    * @throws InvalidAttributeValueException if saltCode not defined
    */
   public void addSaltByCode(OEGraphMol mol, String saltCode, int multiplicity) throws InvalidAttributeValueException {
      assert ! isDeleted;
      
      Salt salt = SaltMap.getSaltCodeMap().getSaltByCode(saltCode);
      if(salt == null)
         throw new InvalidAttributeValueException("Unknown Salt: " + saltCode);
      
      OETools.combineStructures(mol, salt.mol, multiplicity);
   }
   
   /**
    * Add the structure of the salt (saltCode) to the molecule.
    * 
    * @param multiplicity number of copies to add to the molecule, the multiplicity 
    *    in the salt table is ignored.
    *    
    * @throws InvalidAttributeValueException if saltCode not defined
    */
   public String addSaltByCode(String molStr, String saltCode, int multiplicity) throws InvalidAttributeValueException {
      assert ! isDeleted;
      
      mol.Clear();
      molStrToMol(mol, molStr);
      
      addSaltByCode(mol, saltCode, multiplicity);
      return OETools.molToString(mol);
   }

   
  
   /**
    * Add the structure of the salt (saltName) to the molecule.
    * 
    * This will add multiple copies if the salt has a multiplicity.
    * 
    * @throws InvalidAttributeValueException if saltCode not defined
    */
   public void addSaltByName(OEGraphMol mol, String saltName) throws InvalidAttributeValueException {
      assert ! isDeleted;
      
      Salt salt = SaltMap.getSaltCodeMap().getSaltByName(saltName);
      if(salt == null)
         throw new InvalidAttributeValueException("Unknown Salt: " + saltName);
      
      OETools.combineStructures(mol, salt.mol, salt.multiplicity);
   }
   
   /**
    * Add the structure of the salt (saltName) to the molecule.
    * 
    * This will add multiple copies if the salt has a multiplicity.
    * 
    * @throws InvalidAttributeValueException if saltCode not defined
    */
   public String addSaltByName(String molStr, String saltName) throws InvalidAttributeValueException {
      assert ! isDeleted;
      
      mol.Clear();
      molStrToMol(mol, molStr);
      
      addSaltByName(mol, saltName);
      return OETools.molToString(mol);
   }
   
   
   /**
    * Add the structure of the salt (saltName) to the molecule.
    * 
    * @param multiplicity number of copies to add to the molecule, the multiplicity 
    *    in the salt table is ignored.
    *    
    * @throws InvalidAttributeValueException if saltCode not defined
    */
   public void addSaltByName(OEGraphMol mol, String saltName, int multiplicity) throws InvalidAttributeValueException {
      assert ! isDeleted;
      
      Salt salt = SaltMap.getSaltCodeMap().getSaltByName(saltName);
      if(salt == null)
         throw new InvalidAttributeValueException("Unknown Salt: " + saltName);
      
      OETools.combineStructures(mol, salt.mol, multiplicity);
   }
   
   /**
    * Add the structure of the salt (saltName) to the molecule.
    * 
    * @param multiplicity number of copies to add to the molecule, the multiplicity 
    *    in the salt table is ignored.

    * @throws InvalidAttributeValueException if saltCode not defined
    */
   public String addSaltByName(String molStr, String saltName, int multiplicity) throws InvalidAttributeValueException {
      assert ! isDeleted;
      
      mol.Clear();
      molStrToMol(mol, molStr);
      
      addSaltByName(mol, saltName, multiplicity);
      return OETools.molToString(mol);
   }
   
   
   public void close() {
      if( isDeleted ) return;
      
      mol.delete();
      isDeleted = true;
   }
   
   @Override
   protected void finalize() throws Throwable {
      close();
      
      super.finalize();
   }
   
   private static final Pattern EMPTYMolPattern = Pattern.compile(
      "^.*(\\n\\r*|\\r\\n*).*(\\n\\r*|\\r\\n*).*(\\n\\r*|\\r\\n*)  0  0");

   private final void molStrToMol(OEGraphMol mol, String molStr)
   {  if( EMPTYMolPattern.matcher(molStr).find()) 
         OETools.smiToMol(mol, "*");
      else
         OETools.stringToMol(mol, molStr);
   }
/*
   public static void main(String [] args) throws ParseException, JDOMException, IOException {
      long start = System.currentTimeMillis();
      int nStruct = 0;
      
      // create command line Options object
      Options options = new Options();
      Option opt = new Option("sqlFile",true, "sql-xml file");
      opt = new Option("o",true, "output file");
      opt.setRequired(false);
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

      String outFile = cmd.getOptionValue("o");
      
      
      args = cmd.getArgs();
      
      try {
         PrintStream out = System.out;
         if(outFile != null) out = new PrintStream(outFile);
         
         
         
         OEGraphMol mol = new OEGraphMol();
         String molStr = "";
         String saltCode = ""; //
         if(! "1".equals(saltCode) && ! "32".equals(saltCode)) {
            mol.Clear();
            OETools.stringToMol(mol, molStr);
            Salt salt = codeToSaltMap.get(saltCode);
            OETools.combineStructures(mol, salt.mol, salt.multiplicity);
            molStr = OETools.molToString(mol);
         }
         mol.delete();
         
      } catch (Exception e) {
         throw new Error(e);
      }finally {
         System.err.printf("Exported %d structures in %dsec\n",
               nStruct, (System.currentTimeMillis()-start)/1000);
      }
   }

   private static void exitWithHelp(Options options) {
      HelpFormatter formatter = new HelpFormatter();
      formatter.printHelp( INTROText, options );
      System.exit(1);
   }
*/   
 
}

/** Auto-reloading map of saltCodes to Salt Objects */
class SaltMap {
   private static final SQLStatement SALTQuery = 
      SQLStatement.createFromLocal(SaltMap.class, "saltQuery");

   private static final AtomicReference<SaltMap> saltMapRef
                                             = new AtomicReference<SaltMap>();
   
   private final Map<String,Salt> codeToSaltMap = new HashMap<String, Salt>();
   private final Map<String,Salt> nameToSaltMap = new HashMap<String, Salt>();
   
   private final long loadTime = System.currentTimeMillis();

   
   private SaltMap() {
   }
   
   public static SaltMap getSaltCodeMap() {
      SaltMap saltMap = saltMapRef.get();
      if(saltMap == null) 
         saltMap = readSaltMap();
      
      return saltMap;
   }
   
   public Salt getSaltByCode(String saltCode) {
      if( System.currentTimeMillis() - loadTime > SaltAdder.RELAODTime ) {
         readSaltMap();
      }
      
      return saltMapRef.get().codeToSaltMap.get(saltCode);
   }
         
   public Salt getSaltByName(String saltName) {
      if( System.currentTimeMillis() - loadTime > SaltAdder.RELAODTime ) {
         readSaltMap();
      }
      
      return saltMapRef.get().nameToSaltMap.get(saltName);
   }
         
         
   private static SaltMap readSaltMap() {
      SaltMap newSCMap = new SaltMap();
      
      Selecter sel = Selecter.factory( SALTQuery );
      sel.select();
      while(sel.hasNext()) {
         Record rec = sel.next();
         Salt s = new Salt(rec.getStrg(0), rec.getStrg(1), rec.getInt(2), rec.getStrg(3));
         newSCMap.codeToSaltMap.put(s.code, s);
         newSCMap.nameToSaltMap.put(s.name, s);
      }
      saltMapRef.set(newSCMap);
      sel.close();
      
      return newSCMap;
   }
}

class Salt {
   final String code;
   final int multiplicity;
   final String molStr;
   final OEGraphMol mol;
   final String name;
   
   Salt(String code, String name, int multiplicity, String molStr) {
      this.code = code;
      this.name = name;
      this.multiplicity = multiplicity;
      this.molStr = molStr;
      
      this.mol = new OEGraphMol();
      if( ! "1".equals(code) && ! "32".equals(code)) // parent and unknown
      OETools.stringToMol(mol, molStr);
   }
}
