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

/**
 * Ouput one sdf record with a copy of the candidate molecule for each nearest
 * neighbor that matched.
 *
 *   This duplicates each candidate molecule for each match.
 * @author albertgo
 *
 */
public class NNMultiFinderDuplConsumer implements NNMultiFinderConsumerInterface
{  private final oemolothread ofs;
   private final boolean printIds;
   private final String countSimilarityTheshold;

   public NNMultiFinderDuplConsumer(String outFile, boolean printIds, String countSimilarityTheshold)
   {  this.ofs = new oemolothread(outFile);
      this.printIds = printIds;
      this.countSimilarityTheshold = countSimilarityTheshold;
   }

   
   public void consumeResult(OEMolBase mol, TreeSet<Neighbor> nnSet, List<String> referenceIds, int countSim)
   {  if( countSimilarityTheshold != null )
         oechem.OESetSDData(mol, "NNCount_"+countSimilarityTheshold, Integer.toString(countSim));

      if( nnSet.size() > 0)
      {  int count=0;
         for(Neighbor n : nnSet)
         {  oechem.OESetSDData(mol, "NNSim", String.format("%.4f",n.neighBorSim));
            oechem.OESetSDData(mol, "NNIdx", Integer.toString(n.neighBorIdx));
            oechem.OESetSDData(mol, "NNMatchRank", Integer.toString(++count));
            if( printIds )
               oechem.OESetSDData(mol, "NNId", referenceIds.get(n.neighBorIdx));

            oechem.OEWriteMolecule(ofs, mol);
         }
      }else
      {  oechem.OEWriteMolecule(ofs, mol);
      }
   }


   public void close()
   {  ofs.close();
      ofs.delete();
   }
}
