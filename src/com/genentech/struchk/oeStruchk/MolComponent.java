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

import openeye.oechem.OEMolBase;

import com.genentech.oechem.tools.OETools;

/**
 * Holds a component of a salt, mixture or solvated structure.
 * 
 * @author albertgo
 *
 */
class MolComponent implements Comparable<MolComponent> {
   private final int id;
   private final OEMolBase mol;
   private final String canISmi;
   private int occurenceCount = 1;  // start with one since the creation counts as one
   
   MolComponent(int id, OEMolBase mol) {
      this.id = id;
      this.mol = mol;
      this.canISmi = OETools.molToCanSmi(mol, true);
   }

   void delete() {
      mol.delete();
   }

   public int compareTo(MolComponent other) {
      int lenDif = this.mol.NumAtoms() - other.mol.NumAtoms();
      if( lenDif != 0 ) return lenDif;
      
      // ensure that * salt has prefered status and gets removed last
      if("*".equals(canISmi) && ! "*".equals(other.canISmi)) return 1;
      return this.canISmi.compareTo(other.canISmi);
   }

   int getId() {
      return id;
   }

   OEMolBase getMol() {
      return mol;
   }

   String getCanISmi() {
      return canISmi;
   }

   void incrementOccurenceCount() {
      occurenceCount++;
   }
   
   int getOccurenceCount() {
      return occurenceCount;
   }
}

