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

import com.aestel.utility.Message;
import com.aestel.utility.MessageList;
import com.genentech.oechem.tools.Bond;
import com.genentech.struchk.oeStruchk.OEStruchk.StructureFlag;

/**
 * Remove Stereochemistry from double bonds which match given smarts eg. amidine C=N\H.
 *
 * @author albertgo
 *
 */
public class DoubleBondStereoClean extends AbstractStructureCheck{
   private final OESubSearch subSearch;
   private final String smarts;

   public DoubleBondStereoClean(Element tElem) {
      super(tElem);

      smarts = tElem.getTextTrim();
      subSearch = new OESubSearch(smarts);
      if(! subSearch.IsValid())
         throw new Error("Invalid Smarts " + smarts);

      if(! checkExample())
         throw new Error( String.format("Example %s did not have bond stereo to clean.",
                                        getExampleInput()));
   }

   @Override
   public boolean checkStructure(OEGraphMol in, StructureFlag inStereo,
                                 MessageList msgs) {

      OEMatchBaseIter matchIt = subSearch.Match(in);
      while( matchIt.hasNext() ) {
         OEMatchBase match = matchIt.next();

         OEMatchPairBondIter mPairIt = match.GetBonds();
         while( mPairIt.hasNext() ) {
            OEBondBase bd = mPairIt.next().getTarget();
            if(bd.GetOrder() == 2 ) {
               bd.SetBoolData(OEStruchk.STEREOClearTag, true);

               if(bd.HasStereoSpecified()) {
                  msgs.addMessage(new Message(msgTxt, Message.Level.COMMENT, null));

                  Bond.removeDBStereo(bd);
               }
            }
            match.delete();
         }
         mPairIt.delete();
      }
      matchIt.delete();

      return true;
   }

   @Override
   public void delete() { subSearch.delete(); }
}
