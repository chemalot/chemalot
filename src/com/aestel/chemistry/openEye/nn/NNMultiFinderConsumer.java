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
import openeye.oechem.oechem;
import openeye.oechem.oemolothread;

public class NNMultiFinderConsumer implements NNMultiFinderConsumerInterface
{  private final oemolothread ofs;
   private final boolean printIds;
   private final String countSimilarityTheshold;

   /**
    * @param outFile
    * @param printIds
    * @param countSimilarityTheshold count similars >= countSimilar (null == no output)
    */
   public NNMultiFinderConsumer(String outFile, boolean printIds, String countSimilarityTheshold)
   {  this.ofs = new oemolothread(outFile);
      this.printIds = printIds;
      this.countSimilarityTheshold = countSimilarityTheshold;
   }


   public void consumeResult(OEMolBase mol, TreeSet<Neighbor> nnSet, List<String> referenceIds, int countSim)
   {  StringBuilder sb = new StringBuilder(nnSet.size()*9);
      if( nnSet.size() > 0)
      {  for(Neighbor n : nnSet)
         {  sb.append(String.format("%.4f,",n.neighBorSim));
         }
         sb.setLength(sb.length()-1);
      }
      oechem.OESetSDData(mol, "NNSim", sb.toString());

      sb.setLength(0);
      if( nnSet.size() > 0)
      {  for(Neighbor n : nnSet)
         {  sb.append(n.neighBorIdx).append(',');
         }
         sb.setLength(sb.length()-1);
      }
      oechem.OESetSDData(mol, "NNIdx", sb.toString());

      if( printIds )
      {  sb.setLength(0);
         if( nnSet.size() > 0)
         {  for(Neighbor n : nnSet)
            {  sb.append(referenceIds.get(n.neighBorIdx)).append(',');
            }
            sb.setLength(sb.length()-1);
         }
         oechem.OESetSDData(mol, "NNId", sb.toString());
      }
      
      if( countSimilarityTheshold != null )
         oechem.OESetSDData(mol, "NNCount_"+countSimilarityTheshold, Integer.toString(countSim));
      

      oechem.OEWriteMolecule(ofs, mol);
   }


   public void close()
   {  ofs.close();
      ofs.delete();
   }
}
