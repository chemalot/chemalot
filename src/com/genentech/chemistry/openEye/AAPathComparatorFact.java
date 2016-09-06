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

import openeye.oechem.OEGraphMol;
import openeye.oechem.OEMolBase;

import com.aestel.chemistry.openEye.SimComparator;
import com.aestel.chemistry.openEye.SimComparatorFactory;

/**
 * Factory class for AtomAtomPath comparators.
 * This hides the versioning of the comparateors from the users.
 *
 * To use the default version simply use
 * <code>
 *     AAPathComparatorFact(DEFAULT, DEFAULTVersion);
 * </code>
 * @author albertgo
 *
 */
public class AAPathComparatorFact implements SimComparatorFactory<OEMolBase, OEMolBase, SimComparator<OEMolBase>>
{  public static final int DEFAULTVersion = 8;

   public static enum AAPathCompareType
   {
      DEFAULT(7, 0.0),
      INTERMEDIATE(5, 0.0),
      FAST(4, 0.0),

      FUZZY(7,0.0),
      FUZZYINTERMEDIATE(5,0.0),
      FUZZYFAST(4,0.0);

      AAPathCompareType(int pathLen, double minAtSim)
      {  this.pathLen = pathLen;
         this.minAtSim = minAtSim;
      }

      private final int pathLen;
      private final double minAtSim;

      public int getPathLen()
      {
         return pathLen;
      }
      public double getMinAtSim()
      {  return minAtSim;
      }
   }


   private final int pathLen;
   private final double minAtSim;
   private final int version;
   private final AAPathCompareType type;

   public AAPathComparatorFact(AAPathCompareType type, int version)
   {  this.pathLen = type.getPathLen();
      this.minAtSim = type.getMinAtSim();
      this.type = type;
      this.version = version;
   }

   /** returns new objects which should be deleted separately */
   @Override
   public OEMolBase createComparable(OEMolBase in)
   {  return new OEGraphMol(in); // OEMolBase is directly used for comparison
   }

   @Override
   public SimComparator<OEMolBase> createComparator(OEMolBase mol)
   {  if( type.toString().startsWith("FUZZY") )
      {  if( version == 5) return new IAAPathComparator(mol, FuzzyHeadAtomComputer.INSTANCE, pathLen);
         if( version == 6) return new IAAPathComparatorInt(mol, FuzzyHeadAtomComputer.INSTANCE, pathLen);
         if( version == 7) return new IAAPathComparatorFP(mol, FuzzyHeadAtomComputer.INSTANCE, pathLen);
         if( version == 8) return new IAAPathComparatorChar(mol, FuzzyHeadAtomComputer.INSTANCE, pathLen);

         throw new Error("FUZZY only suuported in version 5 6 7 8");
      }
      if( version == 1 )
         return new AAPathComparator(mol, pathLen, minAtSim);

      if( version == 2 )
         return new AAPathComparator2(mol, pathLen, minAtSim);

      if( version == 3 )
         return new AAPathComparator3(mol, pathLen, minAtSim);

      if( version == 4 )
         return new AAPathComparator4(mol, pathLen, minAtSim);

      if( version == 5 ) // uses long slightly faster on linux than int, results are identical
         return new IAAPathComparator(mol,    HeadAtomComputer.INSTANCE, pathLen);

      if( version == 6 ) // uses int
         return new IAAPathComparatorInt(mol, HeadAtomComputer.INSTANCE, pathLen);

      if( version == 7 ) // uses bitvector, slightly faster on linux than long but results differ due to collisions
         return new IAAPathComparatorFP(mol,  HeadAtomComputer.INSTANCE, pathLen);

      if( version == 8 ) // uses Char, slightly faster on linux than long results need validation
         return new IAAPathComparatorChar(mol,  HeadAtomComputer.INSTANCE, pathLen);

      throw new Error("Unsupported version: available versions are 1 2 3 4 5 6 7 8, deafault is: "
                      + DEFAULTVersion);
   }

   @Override
   public void close()
   {  // nothing to do
   }
}
