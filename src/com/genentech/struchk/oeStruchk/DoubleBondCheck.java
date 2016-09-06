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
import com.genentech.struchk.oeStruchk.OEStruchk.StructureFlag;

/**
 * Structure checking rules which checks for double bonds which were drawn with
 * unclear stereochemistry.
 * 
 * @author albertgo
 *
 */
public class DoubleBondCheck extends AbstractStructureCheck{
   
   public DoubleBondCheck(Element tElem) {
      super(tElem);
      
      if(! checkExample())
         throw new Error( String.format("Example %s did not have invalid bond stereo.",
                                        getExampleInput()));
   }
   
   @Override
   public boolean checkStructure(OEGraphMol in, StructureFlag inStereo,
                                 MessageList msgs) {
      boolean hasError = false;
      OEGraphMol molCopy = new OEGraphMol(in);
      
      oechem.OE3DToBondStereo(molCopy);
      
      OEBondBaseIter bIt = in.GetBonds();
      OEBondBaseIter copyIt = molCopy.GetBonds();
      while(  bIt.hasNext() ) {
         OEBondBase bd = bIt.next();
         OEBondBase bdNew = copyIt.next();

         if(bd.GetOrder() != 2) continue;
         
         if(bd.HasStereoSpecified() && ! bd.IsChiral()) {
            hasError = true;
            msgs.addMessage(new Message(String.format(
                  "Bond %s=%s has cis/trans stereo chemistry specified but is symmetric.",
                  Atom.getAtomName(bd.GetBgn()), Atom.getAtomName(bd.GetEnd())), 
                  Message.Level.ERROR,bd));
         }
            
         // check if 3d perception of stereochemistry recognizes that some bonds 
         // do not have stereo chemistry
         // right now there are no known cases where this applies because on input
         // a linear double bond is recognized as not having stereochemistry
         if(bd.HasStereoSpecified() && ! bdNew.HasStereoSpecified() ) {
            hasError = true;
            msgs.addMessage(new Message(String.format("Bond %s=%s has no stereo chemistry.",
                  Atom.getAtomName(bd.GetBgn()), Atom.getAtomName(bd.GetEnd())), 
                  Message.Level.ERROR,bd));
         }
      }
      copyIt.delete();
      bIt.delete();
      molCopy.delete();
      
      return ! hasError;
   }
   
   @Override
   public void delete() { /* nothing to do */ }
}
