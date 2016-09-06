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

import openeye.oechem.OEAtomBaseIter;
import openeye.oechem.OEGraphMol;

import org.jdom.Element;

import com.aestel.utility.Message;
import com.aestel.utility.MessageList;
import com.genentech.oechem.tools.OETools;
import com.genentech.struchk.oeStruchk.OEStruchk.StructureFlag;

/**
 * Structure checking rules which checks if two atoms where drawn very close to
 * each other.
 *
 * @author albertgo
 *
 */
public class TwoDCheck extends AbstractStructureCheck{

   public TwoDCheck(Element tElem) {
      super(tElem);

      if(! checkExample())
         throw new Error( String.format("Example %s did not have 3D coordinates.",
                                        getExampleInput()));
   }

   @Override
   public boolean checkStructure(OEGraphMol in, StructureFlag inStereo,
                                 MessageList msgs) {
      double[] coor = new double[3];
      double maxX, maxY, maxZ, minX, minY, minZ;
      maxX = Double.NEGATIVE_INFINITY;
      maxY = Double.NEGATIVE_INFINITY;
      maxZ = Double.NEGATIVE_INFINITY;
      minX = Double.POSITIVE_INFINITY;
      minY = Double.POSITIVE_INFINITY;
      minZ = Double.POSITIVE_INFINITY;

      OEAtomBaseIter aIt = in.GetAtoms();
      while( aIt.hasNext() ) {
         in.GetCoords(aIt.next(),coor);

         minX = Math.min(coor[0],minX);
         minY = Math.min(coor[1],minY);
         minZ = Math.min(coor[2],minZ);

         maxX = Math.max(coor[0],maxX);
         maxY = Math.max(coor[1],maxY);
         maxZ = Math.max(coor[2],maxZ);
      }
      aIt.delete();

      double minD = Math.min(Math.min(maxX-minX, maxY-minY), maxZ-minZ);
      if( minD > 0.05 ) {
         msgs.addMessage(new Message("Molecule has 3D coordinates, please enter 2D structure.",
                  Message.Level.ERROR, in));
            return false;
      }
      return true;
   }

   @ Override
   public void delete() { /* nothing to do */ }

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
