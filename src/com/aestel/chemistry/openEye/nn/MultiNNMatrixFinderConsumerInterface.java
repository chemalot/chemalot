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

import openeye.oechem.OEMolBase;

public interface MultiNNMatrixFinderConsumerInterface
{
   /**
    * Output result of neighbor search with one compound
    * @param mol query
    * @param baseMolIdx index of query in input
    * @param nnSet set of nearest neighbors to query
    * @param countSimilar count of nearest neighbors above threshold
    */
   public abstract void consumeResult(OEMolBase mol, int baseMolIdx,
                                      TreeSet<Neighbor> nnSet, int countSimilar)
          throws InterruptedException;

   public abstract void close();

}
