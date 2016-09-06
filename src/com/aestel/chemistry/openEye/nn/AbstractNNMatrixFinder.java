/*
   Copyright 2006-2014 Man-Ling Lee & Alberto Gobbi

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

Contact: aestelSW@gmail.com
*/
package com.aestel.chemistry.openEye.nn;

import java.util.ArrayList;
import java.util.List;

import openeye.oechem.*;

import com.aestel.chemistry.openEye.MultiThreadMatrixAlgortihm;
import com.aestel.chemistry.openEye.SimComparator;
import com.aestel.chemistry.openEye.SimComparatorFactory;
/**
 * common methods to do an all by all comparison of records in one input file.
 *
 * @param <T> comparable type
 * @param <Y> comparator type
 */
abstract public class AbstractNNMatrixFinder<T, Y extends SimComparator<T>> implements MultiThreadMatrixAlgortihm
{  protected final SimComparatorFactory<OEMolBase, T, Y> comparableFact;
   protected List<OEMolBase> mols;
   protected List<Y> comparators;

   /**
    * @param compFact for converting the input into a comparable and or comparator.
    * @param inFile file with records to compare.
    */
   public AbstractNNMatrixFinder(SimComparatorFactory<OEMolBase, T, Y> compFact,
                                 String inFile)
   {  comparableFact = compFact;
      mols = readMolecules(inFile);
      comparators = new ArrayList<Y>(2000);
      for( OEMolBase m : mols)
      {  T cmp = compFact.createComparable(m);
         comparators.add(compFact.createComparator(cmp));
         // cmp should be deleted here
      }
   }

   private static List<OEMolBase> readMolecules(String inFile)
   {  List<OEMolBase> myMols = new ArrayList<OEMolBase>(2000);
      OEMolBase mol = new OEGraphMol();
      oemolistream ifs = new oemolistream(inFile);
      int iCounter = 0;

      while(oechem.OEReadMolecule(ifs, mol))
      {  oechem.OESetSDData(mol, "idx", Integer.toString(iCounter));
         myMols.add(mol);
         iCounter++;
         mol = new OEGraphMol();
      }
      System.err.printf("%d molecules read.\n", iCounter);
      ifs.close();
      ifs.delete();
      mol.delete();

      return myMols;
   }

   @Override
   public int getObjectCount()
   {  return comparators.size();
   }

   @Override
   public void close()
   {  for( OEMolBase m: mols)
         m.delete();
      for( Y c: comparators )
         c.close();
   }
}
