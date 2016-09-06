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
import com.genentech.struchk.oeStruchk.OEStruchk.StructureFlag;

/**
 * Structure checking rules which flags an input compound as bad if it contains
 * an unknown atom type.
 * 
 * @author albertgo
 *
 */
public class AtomTypeCheck extends AbstractStructureCheck {
   private final OEStruchk struChecker;
   
   AtomTypeCheck(OEStruchk checker, Element tElem) {
      super(tElem);
      this.struChecker = checker;
   
      if(! checkExample())
         throw new Error( String.format("Example %s was not flagged as having invalid atom.",
                                 getExampleInput()));
   }
   
   @Override
   public boolean checkStructure(OEGraphMol in, StructureFlag inStereo, MessageList msgs) {
      boolean hasError = false;
      
      if( struChecker.isEmptyMol() )
         return true;
      
      OEAtomBaseIter aIt = in.GetAtoms();
      while( aIt.hasNext() ) {
         OEAtomBase at = aIt.next();
         
         if(at.GetAtomicNum() == 0) {
            msgs.addMessage(new Message(
                  String.format("Atom %d has invalid atom type.", at.GetIdx()+1),
                  Message.Level.ERROR, at));
            hasError = true;
         }
         if( ! "".equals(oechem.OEGetAtomComment( at ))) {
            msgs.addMessage(new Message(
                     String.format("Atom %d has alias (%s).", at.GetIdx()+1, oechem.OEGetAtomComment( at )),
                     Message.Level.ERROR, at));
               hasError = true;
         }
      }
      aIt.delete();
      
      return ! hasError;
   }
   
   @Override
   public void delete() { /* nothing to do */ }
}
