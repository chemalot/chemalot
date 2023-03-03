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

/**
 * consumer for NNFinder printing results in smi1<tab>id2<tab>nnSim format.
 *
 * @author albertgo
 *
 */
public class NNFinderVTConsumer implements NNFinderConsumerInterface
{  private final String idTagName;
   private final MultiThreadedPrinter out;

   public NNFinderVTConsumer(String outFile, String idTagName) throws IOException
   {  out = new MultiThreadedPrinter(outFile);
      try
      {  if( idTagName == null )
            out.println("inSmi\trefIdx2\tSim");
         else
            out.println("inIdx1\trefId2\tSim");
      } catch (InterruptedException e)
      {  throw new Error("output queue full imediatly: Should not happen!");
      }

      this.idTagName = idTagName;
   }

   /**
    * countSimilar is ignored for vTab output
    */
   @Override
   public void consumeResult(OEMolBase mol, double maxSim, int nnIdx, String nnId, int countSimilar)
   {  String id1;
      String id2;
      if( idTagName != null )
      {  id1 = oechem.OEGetSDData(mol, idTagName);
         id2 = nnId;     
      }
      else
      {  id1 = OETools.molToCanSmi(mol, true);
         id2 = Integer.toString(nnIdx);
      }

      try
      {  out.println(id1 + '\t' + id2 + '\t' + String.format("%.4f",maxSim));
      } catch (InterruptedException e)
      {  throw new Error("output queue full immediately: Should not happen!");
      }
   }

   @Override
   public void close()
   {  out.close();
   }
}
