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

import java.io.*;
import java.util.TreeSet;

import openeye.oechem.OEMolBase;
import openeye.oechem.oechem;

/**
 * output NN results in tab separated file.
 *
 * This implementation internally stores the similarity matrix and outputs the
 * result in the close method.
 *
 * @author albertgo
 *
 */
public class MultiNNMatrixFinderTabConsumer  implements MultiNNMatrixFinderConsumerInterface
{  private float[][] outMatrix;
   private final String outFile;
   private final String idTag;
   private String[] ids;

   public MultiNNMatrixFinderTabConsumer(String outFile, String idTag)
            throws IOException
   {  this.outFile = outFile;
      this.idTag = idTag;
   }

   /** must be called exactly once before calling consumeResult */
   public void setMatrixSize(int matrixSize)
   {  outMatrix = new float[matrixSize][matrixSize];

      for( int i=0; i< matrixSize; i++)
         for(int j=0; j<matrixSize-1; j++)
            outMatrix[i][j] = 0f;

      for( int i=0; i<matrixSize; i++)
         outMatrix[i][i] = 1f;

      ids = new String[matrixSize];
   }

   @Override
   /**
    * countSimilar is ignored for the tab output
    */
   public void consumeResult(OEMolBase mol, int baseMolIdx, TreeSet<Neighbor> nnSet, int countSimilar)
            throws InterruptedException
   {  if(idTag != null) ids[baseMolIdx] = oechem.OEGetSDData(mol, idTag);

      for(Neighbor nn: nnSet)
      {  outMatrix[baseMolIdx][nn.neighBorIdx] = (float)nn.neighBorSim;
      }
   }

   @Override
   public void close()
   {  try
      {  PrintStream out;
         if(outFile.startsWith("."))
            out = System.out;
         else
            out = new PrintStream(
                     new BufferedOutputStream(new FileOutputStream(new File(outFile))));

         if( idTag != null )
         {  out.print(idTag);
            for( int i=0; i<outMatrix.length; i++)
            {  out.print('\t');
               out.print(ids[i]);
            }
            out.println();
         }

         for( int i=0; i<outMatrix.length; i++)
         {  if( idTag != null )
            {  out.print(ids[i]);
               out.print('\t');
            }

            for(int j=0; j<outMatrix[i].length-1; j++)
               out.printf("%.4f\t",outMatrix[i][j]);

            out.printf("%.4f\n",outMatrix[i][outMatrix[i].length-1]);
         }
         out.close();
      } catch (IOException e)
      {  throw new Error(e);
      }
   }
}
