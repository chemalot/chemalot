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

import com.aestel.utility.MessageList;

import openeye.oechem.OEGraphMol;

public interface StructureKeeperInterface extends StructureCheckInterface {
   /**
    * Get the molecule at the stage of this keeper.
    */
   public OEGraphMol getMolecule();

   public String getKeeperName();
   
   /**
    * Get the transformation messages recorded for this molecule up to the stage
    * of this keeper.
    * 
    * @return a (possibly empty (length == 0)) list of messages.
    */
   public MessageList getStructureMessages();

}
