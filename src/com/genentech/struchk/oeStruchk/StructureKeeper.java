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

import java.util.List;

import openeye.oechem.OEGraphMol;
import openeye.oechem.oechem;

import org.jdom.Element;

import com.aestel.utility.MessageList;
import com.genentech.struchk.oeStruchk.OEStruchk.StructureFlag;
import com.genentech.struchk.oeStruchk.StructureCheckInterface.HydrogenMode;

/**
 * Keep a copy of a molecule with all its messages at the call of the method
 * {@link #checkStructure}.
 *
 * this is not a check as such but a helper class to store a specific state of a
 * molecule during the normalization stages. eg. a salt free molecule.
 *
 * @author albertgo
 */
@SuppressWarnings("unused")
public class StructureKeeper implements StructureKeeperInterface {

   private final String description;
   private final String name;
   private OEGraphMol mol = null;
   private MessageList msgs = null;
   private final String checkName;

   public StructureKeeper(Element elem) {
      description = elem.getChildTextTrim("description");
      name = elem.getAttributeValue("name");
      checkName = elem.getName();
   }

   @Override
   public boolean checkExample() {
      return true;
   }

   @Override
   public boolean checkStructure(OEGraphMol in, StructureFlag inStereo,
                                 MessageList msgs) {
      mol = new OEGraphMol(in); // keep a copy so nobody can change it

      oechem.OESuppressHydrogens(mol,false,false,true);
      oechem.OEPerceiveChiral(mol);

      this.msgs = new MessageList(msgs);

      return true;
   }

   @Override
   public String getDescription() {
      return description;
   }
   @Override
   public String getDescriptionHTML() {
      return getDescription();
   }

   @Override
   public String getExampleInput() {
      return null;
   }

   /**
    * Get the molecule at the stage of this keeper.
    */
   @Override
   public OEGraphMol getMolecule() {
      return mol;
   }

   /**
    * Get the transformation messages recorded for this molecule up to the stage
    * of this keeper.
    *
    * @return a (possibly empty (length == 0)) list of messages.
    */
   @Override
   public MessageList getStructureMessages() {
      return msgs;
   }

   @Override
   public String getKeeperName() {
      return name;
   }

   /**
    * Delete the molecule stored in this keeper.
    *
    * any use of the molecule object retrieved with getMolecule after calling
    * this method will cause errors and possibly crashes!
    */
   @Override
   public void delete() {
      if( mol != null) mol.delete();
      mol = null;
   }

   @Override
   public HydrogenMode getRequiredHydrogenMode() {
      return HydrogenMode.ANY;
   }

   @Override
   public String getCheckName() {
      return checkName;
   }
}
