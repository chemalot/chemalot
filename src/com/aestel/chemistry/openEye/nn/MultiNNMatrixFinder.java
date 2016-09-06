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
import java.util.TreeSet;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorCompletionService;

import openeye.oechem.OEMolBase;

import com.aestel.chemistry.openEye.SimComparator;
import com.aestel.chemistry.openEye.SimComparatorFactory;
/**
 * A nearest Neighbor Finder which compares a full matrix.
 */
public class MultiNNMatrixFinder<T, Y extends SimComparator<T>> extends AbstractNNMatrixFinder<T, Y>
{  private final List<TreeSet<Neighbor>> neighborSets;
   private final int maxNeighbors;
   private final double minSimilarity;
   private final MultiNNMatrixFinderConsumerInterface consumer;
   private final double countSimilarityTheshold;
   private final boolean printAll;


   /**
    *
    * @param c a consumer for the results. the consumeResult method must be thread save
    *        it is closed upon calling {@link #close()}.
    * @param compFact for converting the input into a comparable and or comparator.
    * @param maxNeighbors if 1 only one neighbor is returned, must be > 0
    * @param minSimilarity minimum similarity for consideration, used if maxNeighbors > 1
    * @param printAll if true print input record even if no neighbors are found
    */
   public MultiNNMatrixFinder(String inFile, MultiNNMatrixFinderConsumerInterface c,
            SimComparatorFactory<OEMolBase, T, Y> compFact,
            int maxNeighbors, double minSimilarity, boolean printAll, double countSimilarityTheshold)
   {  super(compFact, inFile);
      this.consumer = c;
      this.maxNeighbors = maxNeighbors;
      this.minSimilarity = minSimilarity;
      this.countSimilarityTheshold = countSimilarityTheshold;
      this.printAll = printAll;

      neighborSets = new ArrayList<TreeSet<Neighbor>>( mols.size() );
   }

   @Override
   public void submitTask(ExecutorCompletionService<Boolean> completionService, int idx)
   {  completionService.submit(new MultiNNMatrixFind(idx));
   }


   /**
    * Compare all reference compounds to the compound whose index is passed to the constructor.
    * @author albertgo
    *
    */
   class MultiNNMatrixFind implements Callable<Boolean>
   {  final Y baseMol;
      final int baseMolIdx;

      MultiNNMatrixFind(int baseMolIdx)
      {  this.baseMol = comparators.get(baseMolIdx);
         this.baseMolIdx = baseMolIdx;
      }

      @Override
      public Boolean call()
      {  int countSimilar = 0;
         double minKnownSim = 10D;
         TreeSet<Neighbor> nnSet = new TreeSet<Neighbor>();

         // NNSearch
         for( int i=0; i< comparators.size(); i++)
         {  if( i == baseMolIdx ) continue;

            Y c = comparators.get(i);
            double sim = baseMol.similarity(c);

            if( sim >= countSimilarityTheshold ) countSimilar++;
            if( sim < minSimilarity ) continue;

            if( nnSet.size() < maxNeighbors )
            {  if( sim < minKnownSim) minKnownSim=sim;
               nnSet.add(new Neighbor(i, sim));
            } else if( sim > minKnownSim )
            {  nnSet.remove(nnSet.last());
               nnSet.add(new Neighbor(i, sim));
               minKnownSim=nnSet.last().neighBorSim;
            }
         }

         neighborSets.add(nnSet);

         OEMolBase mol = mols.get(baseMolIdx);
         try
         {  if( printAll || nnSet.size() > 0 )
               consumer.consumeResult(mol, baseMolIdx, nnSet, countSimilar);
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
