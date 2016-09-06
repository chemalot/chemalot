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

import java.util.TreeSet;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorCompletionService;

import openeye.oechem.*;

import com.aestel.chemistry.openEye.SimComparator;
import com.aestel.chemistry.openEye.SimComparatorFactory;
/**
 * Finds multiple nearest Neighbors in reference file for each compound in the inFile.
 *
 * @param <T> comparable
 * @param <Y> comparator
 */
public class MultiNNFinder<T, Y extends SimComparator<T>> extends AbstractNNFinder<T, Y>
{  final int maxNeighbors;
   final double minSimilarity;
   final oemolithread ifs;
   final NNMultiFinderConsumerInterface consumer;
   private final double countSimilarityTheshold;
   private final boolean printAll;

   /**
    * Find multiple NN in each compound in inFile.
    * @param consumer a consumer for the results. the consumeResult method must be thread save
    *        it is closed upon calling {@link #close()}.
    * @param compFact for converting the input into a comparable and or comparator.
    * @param refFile file with reference molecules
    * @param inFile file with input molecules each input molecule is copared to
    * all references
    * @param maxNeighbors do not report more neighbors
    * @param minSimilarity do not report if similarity is lower
    * @param printAll if true print input record even if no neighbors are found
    * @param countSimilarityTheshold count compounds with similarity above this threshold
    */
   public MultiNNFinder(String inFile, NNMultiFinderConsumerInterface consumer,
            SimComparatorFactory<OEMolBase, T, Y> compFact,
            String refFile, String idTagName,
            int maxNeighbors, double minSimilarity, boolean printAll, double countSimilarityTheshold)
   {  super(compFact, refFile, idTagName);

      this.ifs = new oemolithread(inFile);
      this.consumer = consumer;
      this.maxNeighbors = maxNeighbors;
      this.minSimilarity = minSimilarity;
      this.countSimilarityTheshold = countSimilarityTheshold;
      this.printAll = printAll;

      if( maxNeighbors < 1 ) throw new Error("maxneighbors must be > 0");
   }


   @Override
   public void submitTask(ExecutorCompletionService<Boolean> completionService)
   {  completionService.submit(new NNFindMulti());
   }

   /** inner class to perform multi threaded nn finding **/
   class NNFindMulti implements Callable<Boolean>
   {  /** read one molecule and compare to all references */
      @Override
      public Boolean call()
      {  OEMolBase mol = new OEGraphMol();
         if( ! oechem.OEReadMolecule(ifs, mol) )
         {  mol.delete();
            return Boolean.FALSE;
         }

         TreeSet<Neighbor> nnSet = new TreeSet<Neighbor>();
         Y comp = comparableFact.createComparator(comparableFact.createComparable(mol));
         double minKnownSim = 10D;
         int countSim = 0;

         // NN search
         for( int i=0; i< reference.size(); i++)
         {  Y c = reference.get(i);
            double sim = c.similarity(comp);
            if( sim >= countSimilarityTheshold ) countSim++;
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

         try
         {  if( printAll || nnSet.size() > 0 )
               consumer.consumeResult(mol, nnSet, referenceIds, countSim);
         } catch (InterruptedException e)
         {  e.printStackTrace();

            Thread.currentThread().interrupt();
         }finally
         {  mol.delete();
            comp.close();
         }

         return Boolean.TRUE;
      }
   }

   @Override
   public void close()
   {  super.close();
      consumer.close();
      ifs.close();
   }
}
