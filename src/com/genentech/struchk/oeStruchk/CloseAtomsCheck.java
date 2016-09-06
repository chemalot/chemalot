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
import com.genentech.oechem.tools.Atom;
import com.genentech.oechem.tools.OETools;
import com.genentech.struchk.oeStruchk.OEStruchk.StructureFlag;

/**
 * Structure checking rules which checks if two atoms where drawn very close to 
 * each other.
 * 
 * @author albertgo
 *
 */
public class CloseAtomsCheck extends AbstractStructureCheck{
   
   private final double minDistancePercent;

   public CloseAtomsCheck(Element tElem) {
      super(tElem);
      
      minDistancePercent = Double.parseDouble(tElem.getAttributeValue("minDistance"));
      
      if(! checkExample())
         throw new Error( String.format("Example %s did not have close atoms.",
                                        getExampleInput()));
   }
   
   @Override
   public boolean checkStructure(OEGraphMol in, StructureFlag inStereo,
                                 MessageList msgs) {
      boolean hasError = false;

      if(in.NumBonds() == 0) return true;
      
      // calculate average bond distance
      double avrgBond = OETools.getAverageBondLength(in);
      
      if(avrgBond <=0.001) {
         msgs.addMessage(new Message("Molecule does not have 2D coordinates.",
               Message.Level.ERROR, in));
         return false;
      }
      // squared minimum distance;
      double minDistance2 = avrgBond * minDistancePercent/100;
      minDistance2 *= minDistance2;
      
      // get array of all atoms
      OEAtomBase[] atoms = new OEAtomBase[in.GetMaxAtomIdx()];
      OEAtomBaseIter aIt = in.GetAtoms(); 
      while( aIt.hasNext() ) {
         OEAtomBase at = aIt.next();
         atoms[at.GetIdx()] = at;
      }
      aIt.delete();
      
      // check distances
      for(int i=0; i<atoms.length; i++) {
         if(atoms[i] == null) continue;
         
         for(int j=i+1; j<atoms.length; j++) {
            if(atoms[j] == null) continue;
            
            if(oechem.OEGetDistance2(in, atoms[i], atoms[j]) < minDistance2) {
               hasError = true;
               msgs.addMessage(new Message(String.format("Atom %s is to close to atom %s, please clean your drawing.",
                     Atom.getAtomName(atoms[i]), Atom.getAtomName(atoms[j])),
                     Message.Level.ERROR,
                           atoms[i]));
            }
         }
      }
               
      return ! hasError;
   }
   
   @Override
   public void delete() { /* nothing to do */ }
}
