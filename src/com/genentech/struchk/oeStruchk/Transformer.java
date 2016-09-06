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

import openeye.oechem.*;

import org.jdom.Element;

import com.aestel.chemistry.depict.DepictHelper;
import com.aestel.chemistry.depict.ImageType;
import com.aestel.utility.Message;
import com.aestel.utility.MessageList;
import com.genentech.oechem.tools.AtomBoolPropertyFunctor;
import com.genentech.oechem.tools.OETools;
import com.genentech.struchk.oeStruchk.OEStruchk.StructureFlag;


/**
 * Structure checking role which transforms an input molecule according to a
 * SMIRKS.
 *
 * @author albertgo
 *
 */
public class Transformer extends AbstractStructureCheck {
   private final OEUniMolecularRxn transform;
   private final String smirks;
   private final boolean makeHExplicit;
   private final boolean recheckExoCyclicBonds;
   private final AtomBoolPropertyFunctor clearStereoTagFunctor =
                  new AtomBoolPropertyFunctor(OEStruchk.STEREOClearTag, true);

   /**
    * Create a Transformer from the xml element.
    */
   Transformer(Element tElem) {
      super(tElem);

      makeHExplicit = "true".equalsIgnoreCase(tElem.getAttributeValue("HExplicit"));
      recheckExoCyclicBonds
         = "true".equalsIgnoreCase(tElem.getAttributeValue("recheckExoCyclicBonds"));

      smirks = tElem.getTextTrim();
      transform = new OEUniMolecularRxn(smirks);
      if(! transform.IsValid()) throw new Error("Invalid Smirks " + smirks);

      if(! checkExample())
         throw new Error( String.format("Example %s was not transformed by %s",
                                        getExampleInput(), smirks));
   }

   @Override
   public boolean checkStructure(OEGraphMol in, StructureFlag inStereo,
                                 MessageList msgs) {

      if(makeHExplicit)
         oechem.OEAddExplicitHydrogens(in);

      // transform input molecule according to smirks
      boolean transformed = false;
      if(transform.constCall(in))
      {  transformed = true;
         msgs.addMessage(new Message(msgTxt, Message.Level.COMMENT, null));
      }

      oechem.OESuppressHydrogens(in,false,false,true);
      oechem.OEPerceiveChiral(in);

      if( recheckExoCyclicBonds && transformed )
      {  in.SetBoolData(OEStruchk.NONChiralStereoAssignedTag, false);

         OEAtomBaseIter atIt = in.GetAtoms(clearStereoTagFunctor);
         while(atIt.hasNext())
            atIt.next().DeleteData(OEStruchk.STEREOClearTag);
         atIt.delete();
      }

      return true;
   }

   @Override
   public String getDescriptionHTML() {
      StringBuilder sb = new StringBuilder(2000);
      sb.append(getDescription()).append("<br/>");
      sb.append(AbstractStructureCheck.encodeHTML(smirks)).append("<br/>");

      if(getExampleInput() == null) return sb.toString();

      String smi = getExampleInput();
      sb.append("Example: ")
        .append(DepictHelper.AESTELDepictHelper.getPublicSmilesImageElement(150,100,ImageType.PNG,smi))
        .append("&gt;&gt;");

      OEGraphMol mol = new OEGraphMol();
      oechem.OEParseSmiles(mol, smi);

      MessageList msgs = new MessageList();
      checkStructure(mol, null, msgs );
      smi = OETools.molToCanSmi(mol, true);
      mol.delete();
      sb.append(DepictHelper.AESTELDepictHelper.getPublicSmilesImageElement(100,100,ImageType.PNG,smi))
        .append("<br/>");

      return sb.toString();
   }

   @Override
   public void delete() {
      transform.delete();
      clearStereoTagFunctor.delete();
   }
}
