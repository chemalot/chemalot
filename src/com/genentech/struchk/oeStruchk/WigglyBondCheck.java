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

import openeye.oechem.OEAtomBase;
import openeye.oechem.OEBondBase;
import openeye.oechem.OEBondBaseIter;
import openeye.oechem.OEBondStereo;
import openeye.oechem.OEGraphMol;
import openeye.oechem.OEProperty;

import org.jdom.Element;

import com.aestel.utility.Message;
import com.aestel.utility.MessageList;
import com.genentech.oechem.tools.Atom;
import com.genentech.oechem.tools.Bond;
import com.genentech.oechem.tools.OETools;
import com.genentech.struchk.oeStruchk.OEStruchk.StructureFlag;

public class WigglyBondCheck extends AbstractStructureCheck {

   public WigglyBondCheck(Element ruleElement) {
      super(ruleElement);

      if(! checkExample())
         throw new Error( String.format("Example %s was not flagged as having a wiggly bond",
                                 getExampleInput()));
   }


   @ Override
   public boolean checkStructure(OEGraphMol in, StructureFlag inStereo, MessageList msgs) {
      OEBondBaseIter bit = in.GetBonds();
      while(bit.hasNext()) {
         OEBondBase b = bit.next();
         if( b.GetIntData(OEProperty.BondStereo) == OEBondStereo.Wavy )
         {  OEAtomBase at1 = b.GetBgn();
            OEAtomBase at2 = b.GetEnd();

            // ensure we know which exeocyclic atoms have stereo
            if( at1.IsChiral() || at2.IsChiral() ||
                at1.GetBoolData(OEStruchk.NONChiralStereoAtomTag) ||
                at2.GetBoolData(OEStruchk.NONChiralStereoAtomTag)   )
            {  msgs.addMessage(new Message(
                     String.format("Wiggly bonds not allowed on stereo atoms %s %s",
                                   Atom.getAtomName(at1), Atom.getAtomName(at2)),
                     Message.Level.ERROR, in));
               bit.delete();
               return false;
            }

            b.DeleteData(OEProperty.BondStereo);
            removeStereoFromDoubleBonds(at1);
            removeStereoFromDoubleBonds(at2);
         }
      }
      bit.delete();

      return true;
   }


   private void removeStereoFromDoubleBonds(OEAtomBase at)
   {  OEBondBaseIter bit = at.GetBonds();
      while( bit.hasNext())
      {  OEBondBase b = bit.next();
         if( b.GetOrder() == 2 && b.HasStereoSpecified() )
            Bond.removeDBStereo(b);
      }
      bit.delete();
   }


   @ Override
   public void delete()
   {  // nothing to do
   }


   @ Override
   public boolean checkExample() {
      // no example defined
      if(getExampleInput() == null) return true;

      OEGraphMol mol = new OEGraphMol();
      OETools.stringToMol(mol, getExampleInput());

      MessageList msgs = new MessageList();
      checkStructure(mol, null, msgs );
      mol.delete();

      // check that there is something to report
      return msgs.countMessages() > 0;
   }
}
