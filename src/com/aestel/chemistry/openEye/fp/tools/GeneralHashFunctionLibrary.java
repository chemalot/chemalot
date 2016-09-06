/*
 **************************************************************************
 *                                                                        *
 *          General Purpose Hash Function Algorithms Library              *
 *                                                                        *
 * Author: Arash Partow - 2002                                            *
 * URL: http://www.partow.net                                             *
 * URL: http://www.partow.net/programming/hashfunctions/index.html        *
 *                                                                        *
 * Copyright notice:                                                      *
 * Free use of the General Purpose Hash Function Algorithms Library is    *
 * permitted under the guidelines and in accordance with the most current *
 * version of the Common Public License.                                  *
 * http://www.opensource.org/licenses/cpl.php                             *
 *                                                                        *
 **************************************************************************
*/

package com.aestel.chemistry.openEye.fp.tools;


public class GeneralHashFunctionLibrary
{  private static long NOSIGN = 0x7fffffffffffffffL;

   /**
    * @return an int value smaller (excl.) than maxVal from the input value.
    */
   public static int toInt(long value, int maxVal)
   {  return (int)((value & NOSIGN) % maxVal);
   }

   /**
    * A simple hash function from Robert Sedgwicks Algorithms in C book. I've
    * added some simple optimizations to the algorithm in order to speed up its
    * hashing process.
    */
   public static long RSHash(CharSequence str)
   {  int b     = 378551;
      int a     = 63689;
      long hash = 0;

      for(int i = 0; i < str.length(); i++)
         hash = hash * a + str.charAt(i);
         a    = a * b;

      return hash;
   }
   /* End Of RS Hash Function */


   /**
    * A bitwise hash function written by Justin Sobel
    */
   public static long JSHash(CharSequence str)
   {  long hash = 1315423911;

      for(int i = 0; i < str.length(); i++)
         hash ^= ((hash << 5) + str.charAt(i) + (hash >> 2));

      return hash;
   }
   /* End Of JS Hash Function */


   /**
    * This hash algorithm is based on work by Peter J. Weinberger of AT&T Bell
    * Labs. The book Compilers (Principles, Techniques and Tools) by Aho, Sethi
    * and Ulman, recommends the use of hash functions that employ the hashing
    * methodology found in this particular algorithm.
    */
   public static long PJWHash(CharSequence str)
   {  long BitsInUnsignedInt = 4 * 8;
      long ThreeQuarters     = (BitsInUnsignedInt  * 3) / 4;
      long OneEighth         = BitsInUnsignedInt / 8;
      long HighBits          = (long)(0xFFFFFFFF) << (BitsInUnsignedInt - OneEighth);
      long hash              = 0;
      long test              = 0;

      for(int i = 0; i < str.length(); i++)
      {  hash = (hash << OneEighth) + str.charAt(i);

         if((test = hash & HighBits)  != 0)
            hash = (( hash ^ (test >> ThreeQuarters)) & (~HighBits));
      }

      return hash;
   }
   /* End Of  P. J. Weinberger Hash Function */


   /**
    * Similar to the PJW Hash function, but tweaked for 32-bit processors. Its
    * the hash function widely used on most UNIX systems.
    */
   public static long ELFHash(CharSequence str)
   {  long hash = 0;
      long x    = 0;

      for(int i = 0; i < str.length(); i++)
      {  hash = (hash << 4) + str.charAt(i);

         if((x = hash & 0xF0000000L) != 0)
            hash ^= (x >> 24);
         hash &= ~x;
      }

      return hash;
   }
   /* End Of ELF Hash Function */


   /**
    * This hash function comes from Brian Kernighan and Dennis Ritchie's book
    * "The C Programming Language". It is a simple hash function using a strange
    * set of possible seeds which all constitute a pattern of 31....31...31 etc,
    * it seems to be very similar to the DJB hash function.
    */
   public static long BKDRHash(CharSequence str)
   {  long seed = 131; // 31 131 1313 13131 131313 etc..
      long hash = 0;

      for(int i = 0; i < str.length(); i++)
         hash = (hash * seed) + str.charAt(i);

      return hash;
   }
   /* End Of BKDR Hash Function */


   /**
    * This is the algorithm of choice which is used in the open source SDBM
    * project. The hash function seems to have a good over-all distribution for
    * many different data sets. It seems to work well in situations where there
    * is a high variance in the MSBs of the elements in a data set.
    */
   public static long SDBMHash(CharSequence str)
   {  long hash = 0;

      for(int i = 0; i < str.length(); i++)
         hash = str.charAt(i) + (hash << 6) + (hash << 16) - hash;

      return hash;
   }
   /* End Of SDBM Hash Function */


   /**
    * An algorithm produced by Professor Daniel J. Bernstein and shown first to
    * the world on the usenet newsgroup comp.lang.c. It is one of the most
    * efficient hash functions ever published.
    */
   public static long DJBHash(CharSequence str)
   {  long hash = 5381;

      for(int i = 0; i < str.length(); i++)
         hash = ((hash << 5) + hash) + str.charAt(i);

      return hash;
   }
   /* End Of DJB Hash Function */


   /**
    * An algorithm proposed by Donald E. Knuth in The Art Of Computer Programming
    * Volume 3, under the topic of sorting and search chapter 6.4.
    */
   public static long DEKHash(CharSequence str)
   {  long hash = str.length();

      for(int i = 0; i < str.length(); i++)
         hash = ((hash << 5) ^ (hash >> 27)) ^ str.charAt(i);

      return hash;
   }
   /* End Of DEK Hash Function */


   public static long BPHash(CharSequence str)
   {  long hash = 0;

      for(int i = 0; i < str.length(); i++)
         hash = hash << 7 ^ str.charAt(i);

      return hash;
   }
   /* End Of BP Hash Function */


   /**
    * The FNV hash, short for Fowler/Noll/Vo in honor of the creators, is a very
    * powerful algorithm that, not surprisingly, follows the same lines as
    * Bernstein's modified hash with carefully chosen constants. This algorithm
    * has been used in many applications with wonderful results, and for its
    * simplicity, the FNV hash should be one of the first hashes tried in an
    * application. It is also recommended that the FNV website  be visited for
    * useful descriptions of how to modify the algorithm for various uses.
    */
   public static long FNVHash(CharSequence str)
   {  long fnv_prime = 0x811C9DC5;
      long hash = 0;

      for(int i = 0; i < str.length(); i++)
      {  hash *= fnv_prime;
         hash ^= str.charAt(i);
      }

      return hash;
   }
   /* End Of FNV Hash Function */

/**
 *
 */
   /**
    * An algorithm produced by me Arash Partow. I took ideas from all of the
    * above hash functions making a hybrid rotative and additive hash function
    * algorithm. There isn't any real mathematical analysis explaining why one
    * should use this hash function instead of the others described above other
    * than the fact that I tired to resemble the design as close as possible to
    * a simple LFSR. An empirical result which demonstrated the distributive
    * abilities of the hash algorithm was obtained using a hash-table with
    * 100003 buckets, hashing The Project Gutenberg Etext of Webster's
    * Unabridged Dictionary, the longest encountered chain length was 7,
    * the average chain length was 2, the number of empty buckets was 4579.
    * Below is a simple algebraic description of the AP hash function:
    */
   public static long APHash(CharSequence str)
   {  long hash = 0xAAAAAAAA;

      for(int i = 0; i < str.length(); i++)
      {  if ((i & 1) == 0)
            hash ^= ((hash << 7) ^ str.charAt(i) * (hash >> 3));
         else
            hash ^= (~((hash << 11) + str.charAt(i) ^ (hash >> 5)));
      }

      return hash;
   }
   /* End Of AP Hash Function */



   public static void main(String args[])
   {  String key = "abcdefghijklmnopqrstuvwxyz1234567890";

      System.out.println("General Purpose Hash Function Algorithms Test");
      System.out.println("By Arash Partow - 2002\n");
      System.out.println("Key: " + key);
      System.out.println(" 1. RS-Hash Function Value:   " + RSHash(key));
      System.out.println(" 2. JS-Hash Function Value:   " + JSHash(key));
      System.out.println(" 3. PJW-Hash Function Value:  " + PJWHash(key));
      System.out.println(" 4. ELF-Hash Function Value:  " + ELFHash(key));
      System.out.println(" 5. BKDR-Hash Function Value: " + BKDRHash(key));
      System.out.println(" 6. SDBM-Hash Function Value: " + SDBMHash(key));
      System.out.println(" 7. DJB-Hash Function Value:  " + DJBHash(key));
      System.out.println(" 8. DEK-Hash Function Value:  " + DEKHash(key));
      System.out.println(" 9. BP-Hash Function Value:   " + BPHash(key));
      System.out.println(" 9. FNV-Hash Function Value:  " + FNVHash(key));
      System.out.println("10. AP-Hash Function Value:   " + APHash(key));
   }
}
