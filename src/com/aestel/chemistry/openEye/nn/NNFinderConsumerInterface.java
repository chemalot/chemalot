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

import openeye.oechem.OEMolBase;

/**
 * Oututs single most similar compound to query
 * @author albertgo
 *
 */
public interface NNFinderConsumerInterface
{

   /**
    * 
    * @param mol    query
    * @param maxSim similarity of nearest neighbor
    * @param nnIdx  index of nearest neighbor
    * @param nnId   id of nearest neighbor
    * @param countSim number of compounds above similarity Threshold
    */
   public abstract void consumeResult(OEMolBase mol, double maxSim, int nnIdx,
            String nnId, int countSim);

   public abstract void close();

}
