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
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

/**
 * Start a separate thread to print lines deposited to a blocking Queue.
 *
 * @author albertgo
 *
 */
public class MultiThreadedPrinter
{  private static final String ENDofFile = "sentinel";

   final PrintStream out;
   BlockingQueue<String> outQueue = new ArrayBlockingQueue<String>(200);
   private Thread printerThread;

   public MultiThreadedPrinter(String outFile) throws IOException
   {  if(outFile.startsWith("."))
         out = System.out;
      else
         out = new PrintStream(
                  new BufferedOutputStream(new FileOutputStream(new File(outFile))));

      printerThread = new Thread( new Runnable() {
         public void run()
         {  String s;
            try
            {  while( (s=outQueue.take()) != ENDofFile)
               {  out.println(s);
               }
            } catch (InterruptedException e)
            {  e.printStackTrace();
            }
         }
      }, "MultiThreadedPrinter");
      printerThread.setDaemon(true);
      printerThread.start();
   }


   public void println(String line) throws InterruptedException
   {  outQueue.put(line);
   }

   public void printfln(String format, Object ... args) throws InterruptedException
   {  outQueue.put(String.format(format, args));
   }

   public void close()
   {  try
      {  outQueue.put(ENDofFile);
         printerThread.join();
      } catch (InterruptedException e)
      {  e.printStackTrace();
      }  // mark end of output
      out.close();
   }
}
