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

import java.util.*;

import openeye.oechem.*;

import org.jdom.Element;

import com.aestel.utility.Message;
import com.aestel.utility.MessageList;
import com.genentech.oechem.tools.Atom;
import com.genentech.struchk.oeStruchk.OEStruchk.StructureFlag;

/**
 * Structure checking rules that checks for invalid atomic valences.
 *
 * Atoms which have not been specified may have any valence.
 *
 * @author albertgo
 *
 */
public class AtomValenceCheck extends AbstractStructureCheck {

   public static final int ACCTABLEFragNum = oechem.OEGetTag("ACCTABLEFragNum");

   private final OESubSearch[] accebtableFrags;

   /** Contains all atomic number for atoms for which the valence should be checked */
   private final Set<Integer> valenceCheckAtomSet = new HashSet<Integer>();

   /** Contains Strings of the format atomNum-charge-valence off all
    *  valid valences of all atoms for which valence check is defined */
   private final Set<String> validValenceSet = new HashSet<String>();

   AtomValenceCheck(Element tElem) {
      super(tElem);

      accebtableFrags = getAcceptableFragments(tElem);

      getAtomValenceConfig(tElem, valenceCheckAtomSet, validValenceSet);

      if(! checkExample())
         throw new Error( String.format("Example %s was not flagged as having invalid atom valence.",
                                 getExampleInput()));
   }


   @Override
   public boolean checkStructure(OEGraphMol in, StructureFlag inStereo, MessageList msgs) {
      boolean hasError = false;

      flagAcceptableFragments(in);

      OEAtomBaseIter aIt = in.GetAtoms();
      while(  aIt.hasNext() ) {
         OEAtomBase at = aIt.next();

         if(at.HasData(ACCTABLEFragNum))
            continue;   // at is part of an acceptable fragment

         int atomNr = at.GetAtomicNum();
         if(! valenceCheckAtomSet.contains(atomNr))
            continue;   // no checks for this atom

         int valence = at.GetValence();
         int charge = at.GetFormalCharge();

         if(validValenceSet.contains(String.format("%d-%d-%d", atomNr, charge, valence)))
            continue;   // valid known valence

         msgs.addMessage(new Message(String.format("Atom %s has invalid valence=%d charge=%d.",
               Atom.getAtomName(at), valence, charge),
               Message.Level.ERROR, at ));

         hasError = true;
      }
      aIt.delete();

      return ! hasError;
   }


   /** Flag atoms which are part of the acceptable fragments using
    *  SetBoolData(ACCTABLEFragNum).
    */
   private void flagAcceptableFragments(OEGraphMol in) {

      for(int frag=0; frag<accebtableFrags.length; frag++) {
         OESubSearch fragSearch = accebtableFrags[frag];

         OEMatchBaseIter matchIt = fragSearch.Match(in);
         while( matchIt.hasNext() ) {
            OEMatchBase match = matchIt.next();

            OEMatchPairAtomIter mPairIt = match.GetAtoms();
            while( mPairIt.hasNext() ) {
               mPairIt.next().getTarget().SetBoolData(ACCTABLEFragNum, true);
            }
            mPairIt.delete();
            match.delete();
         }
         matchIt.delete();
      }
   }


   private OESubSearch[] getAcceptableFragments(Element tElem) {
      List<OESubSearch> accepbtableFrag = new ArrayList<OESubSearch>();

      Element aFragEle = tElem.getChild("acceptableFragments");
      if(aFragEle == null) return new OESubSearch[0];

      for(Object elemO : aFragEle.getChildren("fragment")) {
         Element fragEl = (Element)elemO;

         String smarts = fragEl.getTextTrim();

         OESubSearch subSearch = new OESubSearch(smarts);
         if(! subSearch.IsValid())
            throw new Error("Invalid Smarts " + smarts);
         accepbtableFrag.add(subSearch);
      }

      return accepbtableFrag.toArray(new OESubSearch[accepbtableFrag.size()]);
   }


   private static Set<Integer> getAtomValenceConfig(Element tElem,
         Set<Integer> valenceCheckAtoms, Set<String> validValences) {
      for(Object atElementO :  tElem.getChildren("atom")) {
         Element atElement = (Element)atElementO;

         String atom = atElement.getAttributeValue("symbol").trim();
         int atomNr = oechem.OEGetAtomicNum(atom);

         valenceCheckAtoms.add(atomNr);

         for(Object valenceO : atElement.getChildren("valence")) {
            Element valenceEle = (Element) valenceO;

            String charge = valenceEle.getAttributeValue("charge");
            if(charge == null || (charge = charge.trim()).length() == 0)
               throw new Error("Invalid empty charge for atom: " + atomNr);

            String valenceStr = valenceEle.getAttributeValue("values");
            if(valenceStr == null || (valenceStr = valenceStr.trim()).length() == 0)
               throw new Error("Invalid empty valence for atom: " + atomNr);

            for(String val : valenceStr.split("[ \t,;]"))
               validValences.add(String.format("%d-%s-%s", atomNr, charge, val));
         }
      }

      return valenceCheckAtoms;
   }

   @Override
   public void delete() {
      for(int frag=0; frag<accebtableFrags.length; frag++) {
         OESubSearch fragSearch = accebtableFrags[frag];
         fragSearch.delete();
      }
   }
}
