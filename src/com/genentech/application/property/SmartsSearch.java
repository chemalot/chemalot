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
package com.genentech.application.property;
import openeye.oechem.*;

/* counts the number of times a smarts string is matched in a molecule */
public class SmartsSearch {
   private SmartsSearch()
   {  // only static methods here
   }

   public static int search (OEGraphMol mol, String Smarts) {
      int count =0;
      OESubSearch ss = new OESubSearch(Smarts);
      OEMatchBaseIter match = ss.Match(mol, true);
      while (match.hasNext()) {
         count++;
         match.next();
      }
      return count;
   }
   public static int searchAliphaticRings (OEGraphMol mol, String Smarts) {
      int count =0;
      OESubSearch ss = new OESubSearch(Smarts);
      OEQMolBase tmp = ss.GetPattern();
      int numAtoms = tmp.NumAtoms();
      OEMatchBaseIter match = ss.Match(mol, true);
      while (match.hasNext()) {
         //count++;
         OEMatchBase matchBase = match.next();
         OEAtomBaseIter atomIter = matchBase.GetTargetAtoms();
         boolean isAliphatic=false;
         boolean uniqueRing = false;
         while (atomIter.hasNext()){
            OEAtomBase atom = atomIter.next();
            // at least one atom should not be part of a smaller ring
            // SAME AS: at least one atom's OEAtomGetSmallestRingSize is same as SmartsQuery.numatoms()
            // this eliminates double counting of bicyclic, tricyclic, and other fused rings. 
            if (oechem.OEAtomGetSmallestRingSize(atom) == numAtoms){
               uniqueRing = true;
            }
            
            if ( ! atom.IsAromatic() ) {// at least one atom is not aromatic
               isAliphatic = true;
            }
         }
         if (isAliphatic && uniqueRing) {
            count++;
         }
      }
      return count;
   }
}
