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

import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorCompletionService;

import openeye.oechem.*;

import com.aestel.chemistry.openEye.SimComparator;
import com.aestel.chemistry.openEye.SimComparatorFactory;
/**
 * @param <T> comparable
 * @param <Y> comparator
 */
public class NNFinder<T, Y extends SimComparator<T>> extends AbstractNNFinder<T, Y>
{  final oemolithread ifs;
   final NNFinderConsumerInterface consumer;
   private final double countSimilarityTheshold;


   /**
    * Find NN for each compound in inFile compared to each reference compound.
    *
    * @param consumer a consumer for the results. the consumeResult method must be thread save
    *        the is closed upon calling {@link #close()}.
    * @param compFact for converting the input into a comparable and or comparator.
    * @param refFile file with reference molecules
    * @param inFile file with input molecules each input molecule is copared to
    * all references
    * @param countSimilarityTheshold count compounds with similarity above this threshold
    */
   public NNFinder(String inFile, NNFinderConsumerInterface consumer,
                   SimComparatorFactory<OEMolBase, T, Y> compFact,
                   String refFile, String idTagName, double countSimilarityTheshold)
   {  super(compFact, refFile, idTagName);

      ifs = new oemolithread(inFile);
      this.consumer = consumer;
      this.countSimilarityTheshold = countSimilarityTheshold;
   }

   @Override
   public void submitTask(ExecutorCompletionService<Boolean> completionService)
   {  completionService.submit(new NNFind());
   }

   /**
    * Inner class for multi threaded comparison.
    * Read one moleucle form ifs and comare to all references.
    * @author albertgo
    *
    */
   class NNFind implements Callable<Boolean>
   {  @Override
      public Boolean call()
      {  OEMolBase mol = new OEGraphMol();
         if( ! oechem.OEReadMolecule(ifs, mol) )
         {  mol.delete();
            return Boolean.FALSE;
         }

         int countSim = 0;
         Y comp = comparableFact.createComparator(comparableFact.createComparable(mol));
         double maxSim = -1;
         int nnIdx = -1;
         // NNSearch
         for( int i=0; i< reference.size(); i++)
         {  Y c = reference.get(i);
            double sim = c.similarity(comp);
            if( sim > maxSim )
            {  maxSim=sim;
               nnIdx = i;
            }

            if( sim >= countSimilarityTheshold ) countSim++;
         }

         String id = null;
         if(referenceIds.size() > nnIdx && nnIdx != -1) id = referenceIds.get(nnIdx);

         consumer.consumeResult(mol, maxSim, nnIdx, id, countSim);

         mol.delete();
         comp.close();

         return Boolean.TRUE;
      }
   }


   @Override
   public void close()
   {  super.close();
      ifs.close();
      ifs.delete();
      consumer.close();
   }
}
