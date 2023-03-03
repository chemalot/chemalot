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

import openeye.oechem.OEMolBase;
import openeye.oechem.oechem;

import com.genentech.oechem.tools.OETools;

public class NNMatrixFinderVTConsumer  implements NNMatrixFinderConsumerInterface
{  private final MultiThreadedPrinter out;
   private final String idTag;


   public NNMatrixFinderVTConsumer(String outFile, String idTag) throws IOException
   {  out = new MultiThreadedPrinter(outFile);
      this.idTag = idTag;
      try
      {  if( idTag == null )
            out.println("inSmi\trefIdx2\tSim");
         else
            out.println("inIdx1\trefIdx2\tSim");
      } catch (InterruptedException e)
      {  throw new Error("output queue full imediatly: Should not happen!");
      }
   }


   /**
    * countSimilar is ignored for vTab output
    */
   public void consumeResult(OEMolBase mol, double maxSim, int nnIdx, int countSimilar)
      throws InterruptedException
   {  String id1 = idTag == null ? 
                       OETools.molToCanSmi(mol, true)
                     : oechem.OEGetSDData(mol, idTag);
      String id2 = Integer.toString(nnIdx);
      out.println(id1 + '\t' + id2 + '\t' + String.format("%.4f",maxSim));
   }

   public void close()
   {  out.close();
   }
}
