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
       if( at == null )
           return "null";

      return oechem.OEGetAtomicSymbol(at.GetAtomicNum()) + (at.GetIdx()+1);
   }

   /**
    * Get the string representation of the atom in the following format: atom
    * index, atomic symbol and charge.
    */
   public static final String toString(OEAtomBase at) {
       if (at == null)
           return "null";

       StringBuilder atomString = new StringBuilder();
       atomString.append(String.format("%2d", at.GetIdx() + 1))
               .append(String.format("%-1s", oechem.OEGetAtomicSymbol(at.GetAtomicNum())));

       int charge = at.GetFormalCharge();
       if (charge == 0)
           return String.format("%-5s", atomString.toString());

       if (charge < 0)
           atomString.append("-");
       else
           atomString.append("+");

       if (Math.abs(charge) > 1)
           atomString.append(Math.abs(charge));

       return String.format("%-5s", atomString.toString());
   }
 
   /**
    * Get the string representation of the atom in the following format: atom
    * index, atomic symbol and charge, number of implicit and explicit hydrogens and
    * chirality.
    */
   public static final String toString(OEAtomBase at, boolean props) {
       if (at == null)
           return "null";

       if (props == false)
           return toString(at);

       StringBuilder atomString = new StringBuilder(toString(at));

       atomString.append(" HImp=" + at.GetImplicitHCount())
               .append(" HExp=" + (at.GetTotalHCount() - at.GetImplicitHCount()))
               .append(" Chiral=" + (at.IsChiral() ? "T" : "F"));

       return atomString.toString();
   }
   
   /** Find any atom connected to at1 which is not at2 
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
