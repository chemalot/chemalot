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
import java.util.concurrent.ExecutorCompletionService;

import openeye.oechem.*;

import com.aestel.chemistry.openEye.MultiThreadAlgortihm;
import com.aestel.chemistry.openEye.SimComparator;
import com.aestel.chemistry.openEye.SimComparatorFactory;

/**
 * Common methods for all NearNeighbor Finder
 *
 * @param <T> comparable the type of the objects to be compared
 * @param <Y> comparator the type of the references
 */
public abstract class AbstractNNFinder<T, Y extends SimComparator<T>> implements MultiThreadAlgortihm
{  protected final List<Y> reference;
   protected final ArrayList<String> referenceIds;
   protected final SimComparatorFactory<OEMolBase, T, Y> comparableFact;

   /**
    * @param compFact factory to create comparable and comaparator objects.
    * @param refFile reference file, find near neighbors in this list
    * @param idTagName tag name in reference record that is reported as near neighbor ID
    */
   public AbstractNNFinder(SimComparatorFactory<OEMolBase, T, Y> compFact,
                           String refFile, String idTagName)
   {  reference = new ArrayList<Y>(2000);
      referenceIds = new ArrayList<String>(2000);
      comparableFact = compFact;

      readReferenceFile(refFile, idTagName);
   }


   private void readReferenceFile(String refFile, String idTagName)
   {  OEMolBase mol = new OEGraphMol();
      oemolistream ifs = new oemolistream(refFile);
      int iCounter = 0;

      while(oechem.OEReadMolecule(ifs, mol))
      {  if( ! mol.IsValid() )
            throw new Error("Invalid molecule in refrence file: " + (++iCounter));

         oechem.OESetSDData(mol, "CentIdx", Integer.toString(iCounter));
         iCounter++;
         Y c = comparableFact.createComparator(comparableFact.createComparable(mol));
         if( idTagName != null )
         {  String val = oechem.OEGetSDData(mol, idTagName);
            referenceIds.add(val);
            if( val == null || val.length() == 0)
               System.err.println("empty reference id for reference "
                                 + referenceIds.size());
         }
         reference.add(c);
      }
      System.err.printf("%d reference molecules read.\n", iCounter);
      ifs.close();
      mol.delete();
      ifs.delete();
   }

   @Override
   public void close()
   {  comparableFact.close();
   }

   @Override
   public abstract void submitTask(ExecutorCompletionService<Boolean> completionService);
}
