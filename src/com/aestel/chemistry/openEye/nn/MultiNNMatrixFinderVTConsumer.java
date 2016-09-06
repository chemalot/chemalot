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

import java.io.IOException;
import java.util.TreeSet;

import openeye.oechem.OEMolBase;

import com.genentech.oechem.tools.OETools;

/**
 * Output NN result as a vertical table with the format: idx1 <tab> idx2 <tab> sim
 * @author albertgo
 *
 */
public class MultiNNMatrixFinderVTConsumer implements MultiNNMatrixFinderConsumerInterface
{  private MultiThreadedPrinter out;

   public MultiNNMatrixFinderVTConsumer(String outFile) throws IOException
   {  out = new MultiThreadedPrinter(outFile);
      try
      {  out.println("inSmi\trefIdx2\tSim");
      } catch (InterruptedException e)
      {  throw new Error("output queue full imediatly: Should not happen!");
      }
   }


   /* (non-Javadoc)
    * @see MultiNNMatrixFinderConsumerInterface#consumeResult(openeye.oechem.OEMolBase, java.util.TreeSet)
    * countSimilar is ignored for the VTab output
    */
   @Override
   public void consumeResult(OEMolBase mol, int baseMolIdx, TreeSet<Neighbor> nnSet, int countSimilar)
   throws InterruptedException
   {  for(Neighbor n : nnSet)
      {  String smi1 = OETools.molToCanSmi(mol, true);
         String id2 = Integer.toString(n.neighBorIdx);
         out.println(smi1 + '\t' + id2 + '\t' + String.format("%.4f",n.neighBorSim));
      }
   }

   /* (non-Javadoc)
    * @see com.aestel.chemistry.openEye.nn.MultiNNMatrixFinderConsumerInterface#close()
    */
   @Override
   public void close()
   {  out.close();
   }
}
