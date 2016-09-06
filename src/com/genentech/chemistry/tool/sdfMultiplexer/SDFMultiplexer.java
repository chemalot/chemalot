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
import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.atomic.AtomicInteger;

import com.aestel.io.IOUtil;

public class SDFMultiplexer
{  public static final String ENDOfQueue = "ENDOfQueue";

   private final ArrayList<Worker> workers = new ArrayList<Worker>(20);
   private final BlockingQueue<String> outQue;
   private final BlockingQueue<String> inQue;
   private final BufferedReader in;
   private final Thread outputWriter;
   private final boolean stopOnError;
   private final boolean groupByTitle;
   private final boolean groupByAtomCount;
   private final long nBatch;
   private final int nMaxCrashRecover;
   private final long startDelay;
   private final long startTime;
   private final AtomicInteger nCrashes;


   public SDFMultiplexer(BufferedReader in, PrintWriter out, int nMaxWorker,
                         long nBatch, int nCrashRevovery, long startDelay,
                         boolean stopOnError, boolean groupByTitle, boolean groupByAtomCount)
   {  inQue = new ArrayBlockingQueue<String>((nMaxWorker*3)/2);
      outQue = new ArrayBlockingQueue<String>(nMaxWorker);
      this.in = in;
      outputWriter = new OutputWriterThread(outQue, out);
      this.nBatch = nBatch;
      this.nMaxCrashRecover = nCrashRevovery;
      this.stopOnError = stopOnError;
      this.groupByTitle = groupByTitle;
      this.groupByAtomCount = groupByAtomCount;
      this.startTime = System.currentTimeMillis();
      this.startDelay = startDelay;
      this.nCrashes = new AtomicInteger(0);
   }

   /**
    * Add nProcesses execution of commandLine to this Multiplexer.
    *
    * This can be called multiple times in case the command line is to be executed
    * on different hosts and eg. the first part of cmdLine contains a ssh command
    * including the hostname.
    *
    * Input and output from cmdLine must be in sdf records.
    *
    * @param cmdLine a command line which must take input from stdin and write
    *                output to stdout.
    * @param errFilePrefix
    */
   private void addWorkers(String shell, String cmdLine,
                           String errFilePrefix, int nProcesses) throws IOException
   {  for(int i=0; i<nProcesses; i++)
      {  Worker w = new Worker(this, workers.size(), shell, cmdLine,
                               groupByTitle, groupByAtomCount, stopOnError,
                               errFilePrefix, nBatch, nMaxCrashRecover);
         workers.add(w);
      }
   }


   private void run()
   {  outputWriter.start();

      for( Worker w : workers)
      {  w.start();

         try // delay execution of workers in case there is a millisec dependence
         {  if( startDelay > 0 )
               Thread.sleep(startDelay);
         } catch (InterruptedException e)
         {  throw new Error(String.format("W%d: Interrupt while sleeping on thread: %s\n",
                               w.getWorkerNum(), e.getMessage()));
         }
      }


      int molCount = 0;
      String line;
      StringBuilder sb = new StringBuilder(4000);

      try
      {  Thread.sleep(100); // allow workers to startup and start reading and
                            // there is not a single worker getting all the fish
         int lineNum=0;
         String lastTitle = "";
         String lastAtBonds = "";
         String title = "";
         String atBonds = "";
         int lastRecordStart = 0;

         while ( (line = in.readLine()) != null )
         {  lineNum++;

            if( lineNum == 1 )
            {  lastRecordStart = sb.length();
               title = line;
            }

            if(groupByAtomCount && lineNum == 4)
               atBonds = line.substring(0,Math.min(6, line.length()));

            sb.append(line).append('\n');

            if( line.equals("$$$$") )
            {  if( groupByTitle )
               {  if( ! lastTitle.equals(title) )
                  {  // output all but last record
                     if( lastRecordStart > 0 )
                     {  inQue.put(sb.substring(0, lastRecordStart));
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
                     {  inQue.put(sb.substring(0, lastRecordStart));
                        String lastRecord = sb.substring(lastRecordStart);
                        sb.setLength(0);
                        sb.append(lastRecord);
                     }

                     lastTitle = title;
                     lastAtBonds = atBonds;
                  }

               }else // not groupByTitle not groupByAtomCount
               {  inQue.put(sb.toString());
                  sb.setLength(0);
               }

               lineNum = 0;
               molCount++;
            }
         }
         if( ( groupByTitle || groupByAtomCount )
                  && sb.length() != 0 && lineNum == 0 )
         {  inQue.put(sb.toString());  // put last group to queue

         }else if( lineNum != 0 )
            System.err.printf("Incomplete molfile at end of input:%s\n", sb.toString());

         in.close();

      } catch (IOException e)
      {  System.err.printf("Error reading input: %s\n", e.getMessage());

      } catch (InterruptedException e)
      {  System.err.printf("Unexpected interrupt (%s) while queueing input mol:%s\n",
                           e.getMessage(), sb.toString());
      }

      try
      {  for(int i=0; i<workers.size(); i++)
            inQue.put(ENDOfQueue);

      } catch (InterruptedException e)
      {  System.err.printf("Unexpected interrupt (%s) while queueing ENDofQueue.\n",
                           e.getMessage());
      }

      for(Worker w: workers)
      {  try
         {  w.join();

         } catch (InterruptedException e)
         {  System.err.printf("W%d: unexpected interrupt (%s) while waiting for worker to complete.\n",
                  w.getWorkerNum(), e.getMessage());
         }
      }

      try
      {  outQue.put(ENDOfQueue);
         outputWriter.join();
      } catch (InterruptedException e)
      {  System.err.printf("Unexpected interrupt (%s) while waiting for the outputWriter.\n",
               e.getMessage());
      }

      System.err.printf("sdfMultiplexer completed. %d molecules read in %dsec (%d crashRecoveries)\n",
                        molCount, (System.currentTimeMillis() - startTime)/1000, nCrashes.intValue());
   }


   /** return queue which holds the molecules to be worked on*/
   public BlockingQueue<String> getInQueue()
   {  return inQue;
   }

   /** return queue which holds the completed molecules*/
   public BlockingQueue<String> getOutQueue()
   {  return outQue;
   }

   void incrementCrashCount()
   {  nCrashes.incrementAndGet();
   }

   public static void main(String...args) throws IOException
   {  ArrayList<String> shells   = new ArrayList<String>(args.length);
      ArrayList<Integer> nProcs = new ArrayList<Integer>(args.length);

      String cmdLine = "";
      String inFName = "";
      String outFName = "";
      String errFilePrefix = "";
      long   startDelay = 0;

      long nBatch = Long.MAX_VALUE;
      int nCrashRevovery = 10;

      String shell = System.getenv("SHELL");
      if( shell != null ) shell = shell + " -fc";

      boolean stopOnError = false;
      boolean groupByAtomCount = false;
      boolean groupByTitle = false;
      int totalProc = 0;
      for(int i=0; i < args.length; i++)
      {  if("-nProc".equals(args[i]))
         {  if(i+1 == args.length)
               exitWithError("-nProc needs one parameters:\n");

            if( shell == null )
               exitWithError("No SHELL enviroment variable, -shell must apear before -nProc.");

            int nProc = Integer.parseInt(args[++i]);
            totalProc += nProc;
            nProcs.add(nProc);
            shells.add(shell);

         }else if( "-shell".equals(args[i]))
         {  if(i+1 == args.length)
               exitWithError("-shell needs one parameters:\n");

            shell = args[++i];

         }else if( "-cmd".equals(args[i]))
         {  if(i+1 == args.length)
               exitWithError("-cmd needs one parameters:\n");

            if( cmdLine.length() > 0 )
               exitWithError("You may not specify multiple: -cmd or -cmdFile");

            cmdLine = args[++i];

         }else if( "-cmdFile".equals(args[i]))
         {  if( cmdLine.length() > 0 )
               exitWithError("You may not specify multiple: -cmd or -cmdFile");

            if(i+1 == args.length)
               exitWithError("-cmdFile needs one parameters:\n");

            cmdLine = IOUtil.fileToString(args[++i]);

         }else if( "-in".equals(args[i]))
         {  if(i+1 == args.length)
               exitWithError("-in needs one parameters:\n");

            if( inFName.length() > 0 )
               exitWithError("-in only one input file allowed:\n");

            inFName = args[++i];

         }else if( "-out".equals(args[i]))
         {  if(i+1 == args.length)
               exitWithError("-out needs one parameters:\n");

            if( outFName.length() > 0 )
               exitWithError("-out only one output file allowed:\n");

            outFName = args[++i];

         }else if( "-batchSize".equals(args[i]))
         {  if(i+1 == args.length)
               exitWithError("-batchSize needs batchCount:\n");

            nBatch = Long.parseLong(args[++i]);

         }else if( "-maxCrashRecover".equals(args[i]))
         {  if(i+1 == args.length)
               exitWithError("-maxCrashRecover needs count:\n");

            nCrashRevovery = Integer.parseInt(args[++i]);

         }else if( "-errFilePrefix".equals(args[i]))
         {  if(i+1 == args.length)
               exitWithError("-errFilePrefix needs prefix:\n");

            errFilePrefix = args[++i];

         }else if( "-stopOnError".equals(args[i]))
         {  stopOnError = true;

         }else if( "-groupByAtomCount".equals(args[i]))
         {  groupByAtomCount = true;

         }else if( "-groupByTitle".equals(args[i]))
         {  groupByTitle = true;

         }else if( "-startDelay".equals(args[i]))
         {  if(i+1 == args.length)
               exitWithError("-startDelay needs int values:\n");

            startDelay = Long.parseLong(args[++i]);

         }else
         {  exitWithError("Unknown option: %s\n", args[i]);
         }
      }
      StringBuilder err = new StringBuilder();
      if( outFName.length() == 0 )
         err.append("-out is required.\n");

      if( inFName.length() == 0 )
         err.append("-in is required.\n");

      if( totalProc <= 0)
         err.append("-nProc is required.\n");

      if( cmdLine.length() == 0 )
         err.append("-cmdLine is required.\n");

      if( err.length() > 0 )
         exitWithError(err.toString());

      BufferedReader in;
      if( inFName.equalsIgnoreCase(".sdf"))
         in = new BufferedReader(new InputStreamReader(System.in));
      else
         in = new BufferedReader(new FileReader(inFName));

      PrintWriter out;
      if( ".sdf".equalsIgnoreCase(outFName))
         out = new PrintWriter(System.out,false);
      else
         out = new PrintWriter(new BufferedWriter(new FileWriter(outFName)));

      SDFMultiplexer multi = new SDFMultiplexer(in, out, totalProc,
                                                nBatch, nCrashRevovery, startDelay,
                                                stopOnError, groupByTitle, groupByAtomCount);

      for(int i=0; i<shells.size(); i++)
        multi.addWorkers(shells.get(i), cmdLine, errFilePrefix, nProcs.get(i));

      multi.run();
   }

   private static void exitWithError(String msgs, String...args)
   {  System.err.printf(msgs, (Object[])args);
      System.err.println(IOUtil.getResource(SDFMultiplexer.class, "SDFMultiplexer.txt"));
      System.exit(1);
   }



   static class OutputWriterThread extends Thread
   {  private final BlockingQueue<String> queue;
      private final PrintWriter out;

      OutputWriterThread(BlockingQueue<String> outQue, PrintWriter out)
      {  this.queue = outQue;
         this.out = out;
      }

      @Override
      public void run()
      {  try
         {  String mol;
            while( (mol=queue.take()) != ENDOfQueue )
            {  out.print(mol);
            }
            out.flush();
            out.close();
         }
         catch (InterruptedException e)
         {  System.err.printf("Unexpected interrupt (%s) while fetching for output.\n",
                  e.getMessage());
         }
      }
   }



}
