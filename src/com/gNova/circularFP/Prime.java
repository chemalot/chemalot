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

package com.gNova.circularFP;



/**
 * Title:        Prime
 * Description:  Prime a static package to compute prime numbers
 * Copyright:    Copyright (c) 2001
 * Company:      None
 * @author       A. Gobbi
 * @version 1.0
 */

public class Prime
{
   private static int[] primes;
   private static int   maxIndex;

   static
   {  primes = new int[50];
      primes[0] = 1;
      primes[1] = 2;
      primes[2] = 3;

      maxIndex = 2;
   }
   /**
    * get the index'th prime number.
    * <br/>
    * getPrime(0) = 1<br/>
    * getPrime(1) = 2<br/>
    * getPrime(2) = 3<br/>
    * ...
    */
   public static int getPrime(int index)
   {  if(maxIndex < index)
      {  computeNthPrime(index);
      }

      return primes[index];
   }

   private synchronized static void computeNthPrime(int index)
   {  int k;
      int newMaxIndex;
      int nextPrime;
      int newPrimes[];

      if(index <= maxIndex) return;    // already computed

      if(index >= primes.length)       // provide more space
      {  newPrimes = new int[index+51];
         System.arraycopy(primes,0,newPrimes,0,maxIndex+1);
      }else
      {  newPrimes = primes;           // fits into current array
      }

      newMaxIndex = maxIndex;
      nextPrime = newPrimes[maxIndex] + 2;   // next odd number
      k = 2;
      do
      {  if(nextPrime % newPrimes[k] == 0)
         {  nextPrime += 2;                     // found divisor
            k = 2;                              // try next odd
            continue;
         }

         if(nextPrime / newPrimes[k] <= newPrimes[k])   // Knuth 1 p. 147
         {  newPrimes[++newMaxIndex] = nextPrime;      // found prime
            if(newMaxIndex == index)
            {  // we are finished and have found our index'th prime
               primes = newPrimes;              // store new array
               maxIndex = newMaxIndex;
               return;
            }

            nextPrime += 2;   // next odd as candidate, find next
            k = 2;
            continue;
         }

         k++;
      }while(true);
   }
   
   public static int getPrimeLargerEqThan(int i)
   {  int idx = 0;
      // we could do a binary search if performance becomes an issue
      while( getPrime(idx) < i)
         idx++;
      
      return getPrime(idx);
   }

   public static void main(String[] args)
   {  System.err.println(getPrime(133));

      for(int i=0; i<400; i++)
         System.err.println(i + "\t" + getPrime(i));
   }
}
