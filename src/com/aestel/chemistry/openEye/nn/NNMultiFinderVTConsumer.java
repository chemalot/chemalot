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
import java.util.List;
import java.util.TreeSet;

import openeye.oechem.OEMolBase;

import com.genentech.oechem.tools.OETools;

/**
 * Output NN result as a vertical table with the format: smi1 <tab> id2 <tab> sim
 * @author albertgo
 *
 */
public class NNMultiFinderVTConsumer implements NNMultiFinderConsumerInterface
{  private final MultiThreadedPrinter out;
   private boolean printIds;

   public NNMultiFinderVTConsumer(String outFile, boolean printIDs) throws IOException
   {  out = new MultiThreadedPrinter(outFile);
      this.printIds = printIDs;

      try
      {  out.println("inSmi\trefIdx2\tSim");
      } catch (InterruptedException e)
      {  throw new Error("output queue full imediatly: Should not happen!");
      }
   }


   /**
    * CountSim is ignored for VTab output
    */
   @Override
   public void consumeResult(OEMolBase mol, TreeSet<Neighbor> nnSet, List<String> referenceIds, int countSim)
      throws InterruptedException
   {  for(Neighbor n : nnSet)
      {  String smi1 = OETools.molToCanSmi(mol, true);
         String id2;
         if( ! printIds )
            id2  = Integer.toString(n.neighBorIdx);
         else
            id2 = referenceIds.get(n.neighBorIdx);

         out.println(smi1 + '\t' + id2 + '\t' + String.format("%.4f",n.neighBorSim));
      }
   }

   @Override
   public void close()
   {  out.close();
   }
}
