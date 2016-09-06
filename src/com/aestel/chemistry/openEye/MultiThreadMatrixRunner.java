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
package com.aestel.chemistry.openEye;

import java.util.concurrent.*;

public class MultiThreadMatrixRunner<T, Y extends SimComparator<T>>
{  private final ExecutorService executor;
   private final ExecutorCompletionService<Boolean> completionService;
   private final int nCpu;
   private final MultiThreadMatrixAlgortihm algorithm;

   /**
    * Run the algorithm on the square matrix of the objects in algorithm by calling
    * the {@link MultiThreadMatrixAlgortihm#submitTask(ExecutorCompletionService, int)}
    * for each column of the matrix.
    *
    * @param alg algorithm to execute, a call to {@link #clone()} will close the algorithm.
    */
   public MultiThreadMatrixRunner(MultiThreadMatrixAlgortihm alg, int nCpu)
   {  this.algorithm = alg;
      this.nCpu = Math.min(nCpu, alg.getObjectCount());
      if( this.nCpu < nCpu )
         System.err.printf("Reduced nCpu to %d, due to lack of input mols.\n", this.nCpu);

      this.executor = Executors.newFixedThreadPool(this.nCpu);
      this.completionService = new ExecutorCompletionService<Boolean>(executor);
   }


   public void run()
   {  long start = System.currentTimeMillis();
      int iCounter = 0;

      // START THREADS
      int nSubmitted = Math.min(nCpu*2, algorithm.getObjectCount());
      for(int i=0; i<nSubmitted; i++)
         algorithm.submitTask(completionService, i);

      // as soon as a threads completes the calculation for one element
      // start the algorithm for the next element
      try
      {  for( int i= nSubmitted; i< algorithm.getObjectCount(); i++)
         {  completionService.take().get();
            algorithm.submitTask(completionService, i);

            iCounter++;
            if(iCounter % (100*nCpu) == 0)
            {  System.err.print(".");
               if(iCounter % (4000*nCpu) == 0)
               {  System.err.printf( " %d %dsec\n",
                        iCounter, (System.currentTimeMillis()-start)/1000);
               }
            }
         }

         // get last nSubmitted to complete
         for(int i=0; i<nSubmitted; i++)
         {  completionService.take().get();
            iCounter++;
         }

      }catch(ExecutionException e)
      {  // call resulted in exception
         // TODO add option to ignore exceptions???
         throw new Error(e);
      }catch(InterruptedException e)
      {  Thread.currentThread().interrupt();
      }

      String alg = algorithm.getClass().getName();
      if( alg.lastIndexOf('.') >= 0) alg = alg.substring(alg.lastIndexOf('.')+1);
      System.err.printf("\n%s: Processed %d structures. nCpu=%d %.1f sec\n", alg,
            iCounter, nCpu, (System.currentTimeMillis()-start)/1000D);
   }

   /**
    * Write molecule to output file.
    *
    * This is thread safe.
    */
   public void close()
   {  algorithm.close();
      executor.shutdown();
   }
}
