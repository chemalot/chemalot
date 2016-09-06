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
package com.genentech.oechem.tools;

import openeye.oechem.*;

/**
 * static final helper methods for OEAtoms
 *
 */
public class Atom {
   private Atom() {}

   public static final void removeChiralInfo(OEAtomBase at) {
      assert at.HasStereoSpecified();
      
      OEAtomBaseVector v = new OEAtomBaseVector();
      OEAtomBaseIter nbriter=at.GetAtoms();
      while( nbriter.hasNext() ) 
         v.add(nbriter.next());
      nbriter.delete();
      
      at.SetStereo(v, OEAtomStereo.Tetrahedral, OEAtomStereo.Undefined);
      v.delete();
   }

   /**
    * Get the atomic name composed of the atomic symbol followed by (index in
    * the molecule+1).
    */
   public static final String getAtomName(OEAtomBase at) {
      return oechem.OEGetAtomicSymbol(at.GetAtomicNum()) + (at.GetIdx()+1);
   }

   /** find any atom connected to at1 which is not at2 
    * @return null if no atom found.
    */
   public static OEAtomBase findExoNeighbor(OEAtomBase at1, OEAtomBase at2) {
      OEAtomBase oldNat = null;
      OEAtomBaseIter nAtIt=at1.GetAtoms(); 
      while( nAtIt.hasNext() ) {
         oldNat = nAtIt.next();
         if( oldNat.GetIdx() != at2.GetIdx() ) break;          // found atom
      }
      nAtIt.delete();
      
      if( oldNat.GetIdx() == at2.GetIdx() ) return null;
      
      return oldNat;
   }
   
}
