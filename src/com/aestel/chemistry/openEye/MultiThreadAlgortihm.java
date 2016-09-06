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

import java.util.concurrent.ExecutorCompletionService;


/**
 * Implement this interface to be able to execute a task in parallel for each
 * record of an SDF file using the {@link MultiThreadRunner}.
 * @author albertgo
 *
 */
public interface MultiThreadAlgortihm
{
   /**
    * Submit a task to the completion service in this method.
    * An example implementation is {@link com.aestel.chemistry.openEye.nn.NNFinder}
    * @param completionService
    */
   public void submitTask(ExecutorCompletionService<Boolean> completionService);

   /**
    * Free any resources.
    */
   void close();

}
