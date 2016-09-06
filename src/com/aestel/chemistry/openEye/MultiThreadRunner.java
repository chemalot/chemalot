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

public class MultiThreadRunner
{  private final ExecutorService executor;
   private final ExecutorCompletionService<Boolean> completionService;
   private final int nCpu;
   private final MultiThreadAlgortihm algorithm;

   /**
    * run the MultiThreadAlgortihm on nCpu's.
    * Execute {@link MultiThreadAlgortihm#submitTask(ExecutorCompletionService)} on nCpus
    * and start additional execution as long as the method returns true.
    *
    * Each call to {@link MultiThreadAlgortihm#submitTask(ExecutorCompletionService)}
    * should calculate one molecule from an input reader.
    *
    * @param alg algorithm to execute, a call to {@link #close()} will close the algorithm.
    */
   public MultiThreadRunner(MultiThreadAlgortihm alg, int nCpu)
   {  this.algorithm = alg;
      this.nCpu = nCpu;
      this.executor = Executors.newFixedThreadPool(nCpu);
      this.completionService = new ExecutorCompletionService<Boolean>(executor);
   }

   public void run()
   {  long start = System.currentTimeMillis();
      int iCounter = 0;

      // START THREADS
      int nSubmitted = nCpu;//+5;
      for(int i=0; i<nSubmitted; i++)
         algorithm.submitTask(completionService);

      // no wait for any task to return and if it returns TRUE resubmit so
      // that it can do more work
      // the task will return false when it finds no more inputs
      while(nSubmitted > 0)
      {  try
         {  Boolean res = completionService.take().get();
            if( res == Boolean.TRUE )
            {  algorithm.submitTask(completionService);

               iCounter++;
               if(iCounter % (100*nCpu) == 0)
               {  System.err.print(".");
                  if(iCounter % (4000*nCpu) == 0)
                  {  System.err.printf( " %d %dsec\n",
                           iCounter, (System.currentTimeMillis()-start)/1000);
                  }
               }
            }else
            {  nSubmitted--;
            }
         }catch(ExecutionException e)
         {  // call resulted in exception
            // TODO add option to ignore exceptions???
            throw new Error(e);
         }catch(InterruptedException e)
         {  Thread.currentThread().interrupt();
         }
      }

      String alg = algorithm.getClass().getName();
      if( alg.lastIndexOf('.') >= 0) alg = alg.substring(alg.lastIndexOf('.')+1);
      System.err.printf("\n%s: Read %d structures. nCpu=%d %.1f sec\n", alg,
            iCounter, nCpu, (System.currentTimeMillis()-start)/1000D);
   }

   public void close()
   {  algorithm.close();
      executor.shutdown();
   }
}
