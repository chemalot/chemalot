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
import openeye.oechem.OESubSearch;

import org.jdom.Element;

import com.aestel.chemistry.depict.DepictHelper;
import com.aestel.chemistry.depict.ImageType;
import com.aestel.utility.Message;
import com.aestel.utility.MessageList;
import com.genentech.struchk.oeStruchk.OEStruchk.StructureFlag;

/**
 * Structure checking rules which flags an input compound as bad if it matches a
 * SMARTS.
 *
 * @author albertgo
 *
 */
public class BadSubstructureCheck extends AbstractStructureCheck {
   private final OESubSearch subSearch;
   private final String smarts;

   BadSubstructureCheck(Element tElem) {
      super(tElem);
      smarts = tElem.getTextTrim();
      subSearch = new OESubSearch(smarts);
      if(! subSearch.IsValid())
         throw new Error("Invalid Smarts " + smarts);

      if(! checkExample())
         throw new Error( String.format("Example %s was not flagged by %s",
                                 getExampleInput(), smarts));
   }

   @Override
   public boolean checkStructure(OEGraphMol in, StructureFlag inStereo,
                                 MessageList msgs) {

      // if smarts of matches this is a bad structure
      if(subSearch.SingleMatch(in))
      {  msgs.addMessage(new Message(msgTxt, Message.Level.ERROR, null));
         return false;
      }
      return true;
   }

   @Override
   public String getDescriptionHTML() {
      StringBuilder sb = new StringBuilder(2000);
      sb.append(getDescription()).append("<br/>");
      sb.append(AbstractStructureCheck.encodeHTML(smarts)).append("<br/>");

      if(getExampleInput() == null) return sb.toString();

      String smi = getExampleInput();
      sb.append("Example: ")
        .append(DepictHelper.AESTELDepictHelper.getPublicSmilesImageElement(150,100,ImageType.PNG,smi))
        .append("<br/>");

      return sb.toString();
   }


   @Override
   public void delete() {
      subSearch.delete();
   }
}
