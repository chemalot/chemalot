/*
   Copyright 2008-2014 Genentech Inc.

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
package com.genentech.chemistry.openEye;

import openeye.oechem.OEExprOpts;
import openeye.oechem.OEGraphMol;
import openeye.oechem.OEMCSMaxAtoms;
import openeye.oechem.OEMCSType;
import openeye.oechem.OEMolBase;

import com.aestel.chemistry.openEye.SimComparator;
import com.aestel.chemistry.openEye.SimComparatorFactory;

public class MCSSComparatorFact implements SimComparatorFactory<OEMolBase, OEMolBase, SimComparator<OEMolBase>>
{  private final MCSSCompareType type;

   public MCSSComparatorFact(MCSSCompareType type)
   {  this.type = type;
   }

   /** returns new objects which should be deleted separately */
   public OEMolBase createComparable(OEMolBase in)
   {  return new OEGraphMol(in); // OEMolBase is directly used for comparison
   }

   public MCSSComparator createComparator(OEMolBase mol)
   {  if( type == MCSSCompareType.DEFAULT )
      {  int type = OEMCSType.Approximate;
         int minAt = 2;
         int atExpr = OEExprOpts.DefaultAtoms;
         int bdExpr = OEExprOpts.DefaultBonds;
         OEMCSMaxAtoms mcsFunc = new OEMCSMaxAtoms();


         return new MCSSComparator(mol, type, atExpr, bdExpr, mcsFunc, minAt);
      }
      throw new Error("Should not happen!");
   }

   public void close()
   {  // nothing to do
   }
}
