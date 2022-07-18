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
 * static final helper methods for OEBonds
 *
 */
public class Bond {
   private Bond() {}

   /**
    * Remove stereo chemistry info from double bond bd.
    * bd must be a double bond having stereo chemistry.
    */
   public static final void removeDBStereo(OEBondBase bd) {
      assert bd.GetOrder() == 2 : "bond is not double bond";
      
      OEAtomBase at1 = bd.GetBgn();
      OEAtomBase at2 = bd.GetEnd();
      OEAtomBase at1N = Atom.findExoNeighbor(at1, at2);
      if(at1N == null) return; // no neighbor found;
      OEAtomBase at2N = Atom.findExoNeighbor(at2, at1);
      if(at2N == null) return; // no neighbor found;
      
      OEAtomBaseVector vec = new OEAtomBaseVector();
      vec.add(at1N);
      vec.add(at2N);
      bd.SetStereo(vec, OEBondStereo.CisTrans, OEBondStereo.Undefined);
      vec.delete();
   }
   
   /**
    * Get the bond representation composed of atomic names and bond type.
    */
   public static final String bondToString(OEBondBase bd)
   {  
       if( bd == null )
           return "null";

       return String.format("%s %s %s",Atom.getAtomName(bd.GetBgn()), 
            bd.GetType()+bd.GetIntType(), Atom.getAtomName(bd.GetEnd()));
   }
   
   /*
   * Get the string representation of a bond (OEBondBase) in the following format:
   * begin atom, bond order symbol, end atom.
   */
   public static final String toString(OEBondBase bd) {
       if( bd == null )
           return "null";

       String atomStrA = Atom.toString(bd.GetBgn(), false);
       String atomStrB = Atom.toString(bd.GetEnd(), false);
       String bondStr = "N/A";

       switch( bd.GetOrder() ) {
           case 1 : bondStr = " - "; break;
           case 2 : bondStr = " = "; break;
           case 3 : bondStr = " # "; break;
       }
       return atomStrA + bondStr + atomStrB;
   }

   /**
    * Get the string representation of the bond in the following format:
    * begin atom, bond order symbol, end atom.
    */   
   public static final String toString(OEBondBase bd, boolean props) { 
       if( bd == null )
           return "null";

       if( props == false )
           return toString(bd);

       StringBuilder bondString = new StringBuilder(toString(bd)); 

       return bondString.toString();
   }

}
