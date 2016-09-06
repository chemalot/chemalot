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
package com.genentech.struchk.oeStruchk;

import openeye.oechem.OEGraphMol;

import com.aestel.utility.MessageList;
import com.genentech.struchk.oeStruchk.OEStruchk.StructureFlag;
/**
 * Interface for structure checking rules.
 * @author albertgo
 */
public interface StructureCheckInterface {

   public enum HydrogenMode {
      SUPRRESSED,
      EXPLICIT,
      ANY
   }
   /**
    * Check the molecule <code>in</code> according to a structure check rule.
    * @param in molecule to be checked, the molecule will potentially be changed.
    * @param inStructFlag structureFlag for molecule, may be null if StruCheckMode = EXTERNAL.
    * @param msgs list to append messages if any.
    * @return true if no Error was found.
    */
   public boolean checkStructure(OEGraphMol in, StructureFlag inStructFlag, MessageList msgs);
   
   public String getDescription();
   public String getDescriptionHTML();
   
   /**
    * Return smiles for an example structure which would fail or be transformed 
    * by this rule.  
    * @return null if a smiles can not represent a rule violation.
    */
   public String getExampleInput();

   /**
    * check the example structure with the rule to verify that it works at least 
    * with the example. 
    */
   public boolean checkExample();
    
    /**
     * returns true if this check needs to have hydrogen suppressed.
     */
   public HydrogenMode getRequiredHydrogenMode();

   /** desructor */ 
   public void delete();
   
   public String getCheckName();
}
