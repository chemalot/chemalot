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

public interface ComponentNormalizerInterface extends StructureCheckInterface {

   /**
    * Remove salts and solvents and duplicate molecule components.
    * 
    * @return false if no components are left, or if multiple components with 
    *         differing molecular formula are left.
    */
   public abstract boolean checkStructure(OEGraphMol in, StructureFlag inStereo,
         MessageList msgs);

   /**
    * return the salt code of the last normalized compound.
    * 
    * @return the salt code or the parent saltCode if the last compound did not 
    *         have counter ions {@link #getParentSaltCode()}.
    */
   public abstract String getSaltCode();

   /**
    * return the salt code of the last normalized compound.
    * 
    * @return the salt code or the parent saltCode if the last compound did not 
    *         have counter ions {@link #getParentSaltCode()}.
    */
   public abstract String getSaltName();
   
   /** return the number of repetitions of the salt molecule in the last Structure */
   public int getSaltCount();
   
   /** return the molecular weight of the salt in the last structure */
   public String getSaltMW();
   
   /** return the molecular formula of the salt in the last structure */
   public String getSaltMF();
   
   /** this returns a singleton object */
   public abstract String getParentSaltCode();

   public abstract void reset();

   public abstract void delete();

}
