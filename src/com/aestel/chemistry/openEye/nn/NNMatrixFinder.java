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

import openeye.oechem.OEMolBase;

import com.aestel.chemistry.openEye.SimComparator;
import com.aestel.chemistry.openEye.SimComparatorFactory;
/**
 * Compute full NxN Matrix of nearest Neighbors for one input set.
 *
 */
public class NNMatrixFinder<T, Y extends SimComparator<T>> extends AbstractNNMatrixFinder<T, Y>
{  final Neighbor[] nearNeighbors;
   final  NNMatrixFinderConsumerInterface consumer;
   private final double countSimilarityTheshold;

   /**
    * Find multiple NN for each compound in inFile compared to others in inFile.
    *
    * @param consumer a consumer for the results. the consumeResult method must be thread save
    *        the is closed upon calling {@link #close()}.
    * @param compFact for converting the input into a comparable and or comparator.
    * @param inFile file with input molecules each input molecule is copared to all others
    * @param countSimilarityTheshold count similars above this threshold
    */
   public NNMatrixFinder(String inFile, NNMatrixFinderConsumerInterface consumer,
            SimComparatorFactory<OEMolBase, T, Y> compFact, double countSimilarityTheshold)
   {  super(compFact, inFile);
      this.consumer = consumer;
      this.countSimilarityTheshold = countSimilarityTheshold;
      nearNeighbors = new Neighbor[comparators.size()];
   }

   @Override
   public void submitTask(ExecutorCompletionService<Boolean> completionService, int idx)
   {  completionService.submit(new NNMatrixFind(idx));
   }


   /** search NN for compound with index passed to constructor
    *
    */
   class NNMatrixFind implements Callable<Boolean>
   {  final Y baseMol;
      final int baseMolIdx;

      NNMatrixFind(int baseMolIdx)
      {  this.baseMol = comparators.get(baseMolIdx);
         this.baseMolIdx = baseMolIdx;
      }

      @Override
      public Boolean call()
      {  double maxSim = -1;
         int nnIdx = -1;
         int countSimilar = 0;

         // NNSearch
         for( int i=0; i< comparators.size(); i++)
         {  if( i == baseMolIdx ) continue;

            Y c = comparators.get(i);
            double sim = baseMol.similarity(c);
            if( sim > maxSim )
            {  maxSim=sim;
               nnIdx = i;
            }

            if( sim >= countSimilarityTheshold )
               countSimilar++;
         }

         if( nnIdx  != -1 )
            nearNeighbors[baseMolIdx] = new Neighbor(nnIdx, maxSim);

         OEMolBase mol = mols.get(baseMolIdx);
         try
         {  consumer.consumeResult(mol, maxSim, nnIdx, countSimilar);
         } catch (InterruptedException e)
         {  e.printStackTrace();
            Thread.currentThread().interrupt();
         }

         return Boolean.TRUE;
      }
   }

   @Override
   public void close()
   {  super.close();
      consumer.close();
   }
}
