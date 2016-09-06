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
public class CleanReactionCenter extends AbstractStructureCheck{

   public CleanReactionCenter(Element tElem) {
      super(tElem);

      if(! checkExample())
         throw new Error( String.format("Example %s did not have Reaction centers. Or centers were not cleaned",
                                        getExampleInput()));
   }

   @Override
   public boolean checkStructure(OEGraphMol in, StructureFlag inStereo,
                                 MessageList msgs) {

      in.SetRxn(false);
      OEAtomBaseIter aIt = in.GetAtoms();
      while( aIt.hasNext() ) {
         OEAtomBase at = aIt.next();
         if( at.GetRxnRole() != OERxnRole.None)
            at.SetRxnRole(OERxnRole.None);
      }
      aIt.delete();

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
      String smi = OETools.molToCanSmi(mol, true);
      mol.delete();

      // check that there is something to report
      return ! smi.contains(">>");
   }
   public static void main(String a[])
   {  OEGraphMol mol = new OEGraphMol();
      OETools.smiToMol(mol,">>CCO");
      System.err.print(OETools.molToString(mol));
      System.err.print(OETools.molToCanSmi(mol, true));

   }
}
