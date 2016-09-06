/*
   Copyright 2008-2015 Genentech Inc.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

*/
package com.genentech.chemistry.tool.sdfMultiplexer;

import java.io.*;
import java.nio.charset.Charset;
import java.util.Arrays;
import java.util.concurrent.BlockingQueue;

public class Worker extends Thread
{
   private static final Charset CHARSet = Charset.forName("US-ASCII");
   private final SDFMultiplexer master;
   final BlockingQueue<String> outQue;
   private final BlockingQueue<String> inQue;
   final int workerNum;
   private final String errFilePrefix;
   private final String[] shell;
   private final String cmdLine;
   private final boolean stopOnError;

   final boolean groupByTitle;
   final boolean groupByAtomCount;
   private long nBatch;
   private int nCrashRecover;

   /**
    * @param master sdfMultiplexer driving this worker.
    * @param workerNum index of this worked in master.
    * @param shell shell to execute command line eg. "" or "ssh hostname"
    * @param cmdLine commandLine to execute $procNum $currentDir will be replaced.
    * @param groupByTitle
    * @param errFilePrefix if not "" write err output to file with this prefix.
    * @param nBatch restart cmdLine after nBatch molecules are completed.
    * @param nCrashRecover maximum number of times to recover from a crash.
    */
   public Worker(SDFMultiplexer master, int workerNum, String shell, String cmdLine,
                 boolean groupByTitle, boolean groupByAtomCount, boolean stopOnError,
                 String errFilePrefix, long nBatch, int nCrashRecover)
   {  this.master = master;
      this.inQue = master.getInQueue();
      this.outQue = master.getOutQueue();
      this.workerNum = workerNum;
      this.nBatch = nBatch;
      this.nCrashRecover = nCrashRecover;
      this.errFilePrefix = errFilePrefix;
      this.groupByTitle = groupByTitle;
      this.groupByAtomCount = groupByAtomCount;
      this.stopOnError = stopOnError;
      this.shell = shell.split("\\s");

      cmdLine = cmdLine.replace("$procNum", Integer.toString(workerNum));
      cmdLine = cmdLine.replace("$currentDir", System.getProperty("user.dir"));
      this.cmdLine = cmdLine;
   }

   /**
    * @param cmdLine
    * @param shell
    * @param workerNum
    * @param batchCount
    * @return
    */
   private String[] compileCommandLine(int workerNum, int batchCount)
   {  String cmdLn = cmdLine.replace("$batchNum", Integer.toString(batchCount));

      String[] cmd = Arrays.copyOf(shell, shell.length+1);
      cmd[cmd.length-1] = cmdLn;
      return cmd;
   }

   @Override
   public void run()
   {  int batchCount = -1;
      boolean completed = false;
      while( !completed && nCrashRecover-- != 0)
      {  batchCount++;
         String[] cmdAndArgs = compileCommandLine(workerNum, batchCount);
         Process process;
         ProgramErrorReaderThread errReaderThread;
         OutputReaderThread outputReaderThread;
         try
         {  StringBuilder sb = new StringBuilder();
            sb.append("starting ").append(workerNum).append(':');
            for(String s: cmdAndArgs) sb.append(" " + s);
            System.err.println(sb.toString());

            process = Runtime.getRuntime().exec(cmdAndArgs);
         } catch (IOException e1)
         {  StringBuilder sb = new StringBuilder();
            for(String s: cmdAndArgs) sb.append(s).append(' ');

            throw new Error(String.format("Error (%s, args) executing commandLine:\n%s\n",
                e1.getMessage(), sb.toString()));
         }

         errReaderThread = new ProgramErrorReaderThread( errFilePrefix, process.getErrorStream() );
         errReaderThread.start();

         outputReaderThread = new OutputReaderThread( process.getInputStream());
         outputReaderThread.start();
         OutputStream procInput = process.getOutputStream();
         String mol = null;
         try
         {  long nBatchLeft = nBatch;

            while( ! completed )
            {  if((mol = inQue.take()) == SDFMultiplexer.ENDOfQueue )
               {  completed = true;
                  break;
               }

               procInput.write(mol.getBytes(CHARSet));

               if( --nBatchLeft == 0)
               {  nCrashRecover++;  // this is not a crash, compensate for decrement
                  nBatchLeft = nBatch;
                  break;
               }
            }

         }catch( IOException e )
         {  // crash
            System.err.printf("W%d: Error (%s) writing mol:%s\n"
                             +"Resubmitting, but some molcules might be missing\n",
                             workerNum, e.getMessage(), mol);
            try
            {  inQue.put(mol);  //requeue failed molecule
            } catch (InterruptedException e1)
            {  System.err.printf("W%d: unexpected interrupt (%s) requeueing mol\n",
                              workerNum, e1.getMessage());
            }

            master.incrementCrashCount();

         } catch (InterruptedException e)
         {  System.err.printf("W%d: unexpected interrupt (%s) queueing mol:%s\n",
                              workerNum, e.getMessage(), mol);

         }finally
         {  try
            {  procInput.close();
            } catch (IOException e)
            {  System.err.printf("W%d: error closing stream to process: %s\n",
                                 workerNum, e.getMessage());
            }

            int exCode = -1;
            try
            {  exCode= process.waitFor();
            } catch (InterruptedException e)
            {  System.err.printf("W%d: unexpected interrupt waiting for process exit: %s\n",
                                 workerNum, e.getMessage());
            }

            try
            {  outputReaderThread.join();
               errReaderThread.join();
            } catch (InterruptedException e)
            {  System.err.printf("W%d: unexpected interrupt waiting for EOF from process: %s\n",
                                 workerNum, e.getMessage());
            }

            System.err.printf("W%d: Shutdown exitCode=%d\n", workerNum, exCode);
            if( exCode != 0 && stopOnError ) nCrashRecover = 0;
         }
      }
      System.err.printf("W%d: finished\n", workerNum);
   }


   public int getWorkerNum()
   {  return workerNum;
   }


   class OutputReaderThread extends Thread
   {  private final BufferedReader reader;

      OutputReaderThread(InputStream outStrm)
      {  reader = new BufferedReader( new InputStreamReader( outStrm ) );
      }


      @Override
      public void run()
      {  String line;
         StringBuilder sb = new StringBuilder(4000);

         try
         {  int lineNum=0;
            String title = "";
            String lastTitle = "";
            String atBonds = "";
            String lastAtBonds = "";
            int lastRecordStart = 0;

            while ( (line = reader.readLine()) != null )
            {  lineNum++;

               if( lineNum == 1 )
               {  lastRecordStart = sb.length();
                  title = line;
               }

               if(groupByAtomCount && lineNum == 4)
                  atBonds = line.substring(0,Math.min(6, line.length()));

               sb.append(line).append('\n');
               if( line.equals("$$$$") )
               {
                  if( groupByTitle )
                  {  if( ! lastTitle.equals(title) )
                     {  // output all but last record
                        if( lastRecordStart > 0 )
                        {  outQue.put(sb.substring(0, lastRecordStart));
                           String lastRecord = sb.substring(lastRecordStart);
                           sb.setLength(0);
                           sb.append(lastRecord);
                        }

                        lastTitle = title;
                        lastAtBonds = atBonds;
                     }

                  } else if( groupByAtomCount )
                  {  if( ! lastAtBonds.equals(atBonds) )
                     {  // output all but last record
                        if( lastRecordStart > 0 )
                        {  outQue.put(sb.substring(0, lastRecordStart));
                           String lastRecord = sb.substring(lastRecordStart);
                           sb.setLength(0);
                           sb.append(lastRecord);
                        }

                        lastTitle = title;
                        lastAtBonds = atBonds;
                     }

                  }else // not groupByTitle not groupByAtomCount
                  {  outQue.put(sb.toString());
                     sb.setLength(0);
                  }

                  lineNum = 0;
               }
            }
            if( ( groupByTitle || groupByAtomCount )
                && sb.length() != 0 && lineNum == 0 )
            {  outQue.put(sb.toString());  // put last group to queue

            }else if( lineNum != 0  )
            {  System.err.printf("W%d: Incomplete molfile at end of output:%s\n",
                                 workerNum, sb.toString());
            }

            reader.close();

         } catch (IOException e)
         {  System.err.printf("w%D: Error reading from output: %s\n", workerNum, e.getMessage());

         } catch (InterruptedException e)
         {  System.err.printf("W%d: unexpected interrupt (%s) while queueing mol:%s\n",
                              workerNum, e.getMessage(), sb.toString());
         }
      }
   }


   class ProgramErrorReaderThread extends Thread
   {  private final BufferedReader reader;
      private final PrintStream err;

      ProgramErrorReaderThread(String errFilePrefix, InputStream errStrm)
      {  reader = new BufferedReader( new InputStreamReader( errStrm ) );


         if( errFilePrefix.length() > 0 )
         {  try
            {  err = new PrintStream(new FileOutputStream(errFilePrefix + '.' + workerNum + ".err.log", true));
            } catch (FileNotFoundException e)
            {  throw new Error("Coold not open erro stream: " + e.getMessage(), e);
            }
         }
         else
         {  err = System.err;
         }
      }


      @Override
      public void run()
      {  String line;
         try
         {  while ( (line = reader.readLine()) != null )
            {  err.println( line );
            }
            reader.close();

         } catch (IOException e)
         {  System.err.printf("W%d: Error reading from STDErr: ", workerNum, e.getMessage());
         }

      }
   }
}
