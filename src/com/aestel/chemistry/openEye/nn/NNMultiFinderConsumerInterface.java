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

import java.util.List;
import java.util.TreeSet;

import openeye.oechem.OEMolBase;

/**
 * output n nearest neighbors to query
 * @author albertgo
 *
 */
public interface NNMultiFinderConsumerInterface
{

   /**
    * 
    * @param mol query
    * @param nnSet set of nearest neighbors
    * @param referenceIds ids of nearest neighbors
    * @param countSim number of compounds with similarity above threshold
    * @throws InterruptedException
    */
   public abstract void consumeResult(OEMolBase mol, TreeSet<Neighbor> nnSet,
                        List<String> referenceIds, int countSim) throws InterruptedException;

   public abstract void close();

}
