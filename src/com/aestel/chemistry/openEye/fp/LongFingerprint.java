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
package com.aestel.chemistry.openEye.fp;

/**
 * This is an implementation of {@link Fingerprint} that stores the bits
 * in a byte array.
 *
 * @author albertgo
 */
final public class LongFingerprint implements Fingerprint
{  private final long[] longs;
   private final int nBits;

   public LongFingerprint(String hex)
   {  // convert from hex to bytes

      long[] l = new long[(hex.length()+15)/16];
      int pos = 0;
      hex = hex.toLowerCase();

      for(int i=0; i<hex.length(); i += 16)
      {  String w = hex.substring(i, Math.min(i+16,hex.length()));
         l[pos++] = myParseLong(w);
      }

      this.longs = l;

      int nb = 0;
      for(long lg : longs)
         nb += Long.bitCount(lg);
      this.nBits = nb;
   }

   private static long myParseLong(String hexStr)
   {  // fix Long.parseLong not being able to deal with rihtmost bit set
      // supposedly fixed in 1.8
      if( hexStr.length() < 16 ) return Long.parseLong(hexStr, 16) << (16 - hexStr.length()) * 4;
      if( hexStr.charAt(0) < '8' )  return Long.parseLong(hexStr, 16);
      hexStr = ((char)Integer.parseInt(hexStr.substring(0,1),16)-8) + hexStr.substring(1,16);
      return Long.parseLong(hexStr, 16) | 0x8000000000000000L;
   }

   public LongFingerprint(long[] longs)
   {  this.longs = longs;
      int nb = 0;
      for(long l : longs)
         nb += Long.bitCount(l);
      this.nBits = nb;
   }

   @Override
   public double tanimoto(Fingerprint other)
   {  if( ! (other instanceof LongFingerprint) )
         throw new Error("Not supported yet needs ByteFingerprint");

      if( this.nBits == 0 && other.getNBits() == 0) return 1;

      long[] fp_1 = this.longs;
      long[] fp_2 = ((LongFingerprint) other).longs;

      int andBitCount = 0;  // number of bits present in both

      int i=0;
      int len;
      len = Math.min(fp_1.length, fp_2.length);

      // we could possibly convert to long in 16 ops and then do the and/or
      // and count on the long values more efficiently but I am not sure
      for( ; i<len; i++)
      {  andBitCount += Long.bitCount(fp_1[i] & fp_2[i]);
      }

      return ( (double) andBitCount )/(this.nBits+other.getNBits()-andBitCount);
   }

   /**
    * modified tanimoto which should be less dependent on the size difference of
    * cf. AtomatomPath Paper
    */
   public double mtanimoto(Fingerprint other)
   {  if( ! (other instanceof LongFingerprint) )
         throw new Error("Not supported yet needs ByteFingerprint");

      if( this.nBits == 0 && other.getNBits() == 0) return 1;

      long[] fp_1 = this.longs;
      long[] fp_2 = ((LongFingerprint) other).longs;

      int andBitCount = 0;

      int len = Math.min(fp_1.length, fp_2.length);

      // we could possibly convert to long in 16 ops and then do the and/or
      // and count on the long values more efficiently but I am not sure
      for(int i=0; i<len; i++)
      {  andBitCount += Long.bitCount(fp_1[i]  & fp_2[i]);
      }


//System.err.printf("and=%d or=%d ", andBitCount, orBitCount);
      return ( (double) andBitCount )/(Math.max(this.nBits, other.getNBits()) * 2 - andBitCount);
   }


   @Override
   public Fingerprint fold(int size)
   {
      throw new Error("fold not implemented yet ByteFingerprint");
   }

   @Override
   public String getBinString()
   {  StringBuilder sb = new StringBuilder(longs.length*64);
      for( long l : longs)
      {  String s = Long.toBinaryString(l);
         for(int i=64-s.length(); i>0; i--) sb.append('0');
         sb.append(s);
      }

      return sb.toString();
   }

   @Override
   public int[] getBits()
   {
      throw new Error("getBits not implemented yet for LongFingerprint");
   }

   @Override
   public String getHexString()
   {  StringBuilder sb = new StringBuilder(longs.length*8);
      for( long l : longs)
         sb.append(String.format("%016x", l));

      return sb.toString();
   }

   @Override
   public int getNBits()
   {  return nBits;
   }


   public static void main(String ... args)
   {  LongFingerprint fp1 = new LongFingerprint("fffffaff8b10979880014b000803110010800040002060081000000002000000000000000c0410100400081e08007000070000000004068020010804004912040010008220180810000000000000008808000000000020400200800000200000000020000000000000000000000000000000000000000000000000000000000000000400800000000002000000000000000050000200000080000000000000000000020000800000000000008000000000000000000000000000000000000000000400000004000000000000000000000000000000000000004004000000000000000000000108008020000000000000000000000021");
      LongFingerprint fp2 = new LongFingerprint("fffffaff2b821f9880214a0088038100100100c004216008104000000200100000000000040010100400081e08007000070000000004068020010804004912040011000220080810000000000002008808040000000020400200800000200000000020008000004000000000000000000000000000000000000000000000000000000400800000000002000000000000000050000200000080000000000000000000020000800000000000008000000000000000000000000000000000000000000400000004000000000000000000000000000000000000004004000000000000000000000108008020000000000000000000000021");

      System.err.printf("n1=%d n2=%d sim=%f\n", fp1.getNBits(), fp2.getNBits(), fp1.tanimoto(fp2));
   }
}
