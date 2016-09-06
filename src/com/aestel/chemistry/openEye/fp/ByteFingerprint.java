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
final public class ByteFingerprint implements Fingerprint
{  private final byte[] bytes;
   private final int nBits;

   public ByteFingerprint(String hex)
   {  // convert from hex to bytes

      byte[] b = new byte[(hex.length()+1)/2];
      int pos = 0;
      hex = hex.toLowerCase();

      for(int i=0; i<hex.length(); i++)
      {  int nible = hex.charAt(i);

         if( nible > '9' )
            nible = nible - ('a' - 10);
         else
            nible = nible - '0';
         assert nible >=0 && nible < 16;

         byte byt = (byte)(nible << 4);
         if( ++i < hex.length() )
         {  nible = hex.charAt(i);
            if( nible > '9' )
               nible = nible - ('a' - 10);
            else
               nible = nible - '0';
            assert nible >=0 && nible < 16;

            byt = (byte)(byt | nible);
         }
         b[pos++] = byt;
      }

      this.bytes = b;

      int nb = 0;
      for(byte bi : bytes)
         nb += BYTEBitcount[bi & 0xFF];
      this.nBits = nb;
   }

   public ByteFingerprint(byte[] bytes)
   {  this.bytes = bytes;
      int nb = 0;
      for(byte b : bytes)
         nb += BYTEBitcount[b & 0xFF];
      this.nBits = nb;
   }

   @Override
   public double tanimoto(Fingerprint other)
   {  if( ! (other instanceof ByteFingerprint) )
         throw new Error("Not supported yet needs ByteFingerprint");

      if( this.nBits == 0 && other.getNBits() == 0) return 1;

      byte[] fp_1 = this.bytes;
      byte[] fp_2 = ((ByteFingerprint) other).bytes;

      int andBitCount = 0;  // number of bits present in both
      int orBitCount = 0;   // number of bits present in either

      int i=0;
      byte[] longerFP;
      int len;
      if( fp_1.length < fp_2.length)
      {  longerFP = fp_2;
         len = fp_1.length;
      }else
      {  longerFP = fp_1;
         len = fp_2.length;
      }

      // we could possibly convert to long in 16 ops and then do the and/or
      // and count on the long values more efficiently but I am not sure
      for( ; i<len; i++)
      {  andBitCount += BYTEBitcount[ fp_1[i]  & fp_2[i]  & 0xFF ];
         orBitCount  += BYTEBitcount[ (fp_1[i] | fp_2[i]) & 0xFF ];
      }

      len = longerFP.length;
      for( ; i<len; i++)
         orBitCount += BYTEBitcount[ longerFP[i] & 0xFF ];

//System.err.printf("and=%d or=%d ", andBitCount, orBitCount);
      return ( (double) andBitCount )/orBitCount;
   }

   /**
    * modified tanimoto which should be less dependent on the size difference of
    * cf. AtomatomPath Paper
    */
   public double mtanimoto(Fingerprint other)
   {  if( ! (other instanceof ByteFingerprint) )
         throw new Error("Not supported yet needs ByteFingerprint");

      if( this.nBits == 0 && other.getNBits() == 0) return 1;

      byte[] fp_1 = this.bytes;
      byte[] fp_2 = ((ByteFingerprint) other).bytes;

      int andBitCount = 0;

      int len = Math.min(fp_1.length, fp_2.length);

      // we could possibly convert to long in 16 ops and then do the and/or
      // and count on the long values more efficiently but I am not sure
      for(int i=0; i<len; i++)
      {  andBitCount += BYTEBitcount[ fp_1[i]  & fp_2[i]  & 0xFF ];
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
   {  StringBuilder sb = new StringBuilder(bytes.length*8);
      for( Byte b : bytes)
      {  String s = Integer.toBinaryString(b & 0xFF);
         for(int i=8-s.length(); i>0; i--) sb.append('0');
         sb.append(s);
      }

      return sb.toString();
   }

   @Override
   public int[] getBits()
   {
      throw new Error("getBits not implemented yet ByteFingerprint");
   }

   @Override
   public String getHexString()
   {  StringBuilder sb = new StringBuilder(bytes.length*8);
      for( Byte b : bytes)
         sb.append(String.format("%02x", b));

      return sb.toString();
   }

   @Override
   public int getNBits()
   {  return nBits;
   }

   /** for each of the 256 byte values store the number of bits set to 1 */
   static final int[] BYTEBitcount = new int[]
   {  0,1,1,2,1,2,2,3, 1,2,2,3,2,3,3,4,  1,2,2,3,2,3,3,4, 2,3,3,4,3,4,4,5,
      1,2,2,3,2,3,3,4, 2,3,3,4,3,4,4,5,  2,3,3,4,3,4,4,5, 3,4,4,5,4,5,5,6,
      1,2,2,3,2,3,3,4, 2,3,3,4,3,4,4,5,  2,3,3,4,3,4,4,5, 3,4,4,5,4,5,5,6,
      2,3,3,4,3,4,4,5, 3,4,4,5,4,5,5,6,  3,4,4,5,4,5,5,6, 4,5,5,6,5,6,6,7,
      1,2,2,3,2,3,3,4, 2,3,3,4,3,4,4,5,  2,3,3,4,3,4,4,5, 3,4,4,5,4,5,5,6,
      2,3,3,4,3,4,4,5, 3,4,4,5,4,5,5,6,  3,4,4,5,4,5,5,6, 4,5,5,6,5,6,6,7,
      2,3,3,4,3,4,4,5, 3,4,4,5,4,5,5,6,  3,4,4,5,4,5,5,6, 4,5,5,6,5,6,6,7,
      3,4,4,5,4,5,5,6, 4,5,5,6,5,6,6,7,  4,5,5,6,5,6,6,7, 5,6,6,7,6,7,7,8
   };

   public static void main(String ... args)
   {  ByteFingerprint fp1 = new ByteFingerprint("fffffaff8b10979880014b000803110010800040002060081000000002000000000000000c0410100400081e08007000070000000004068020010804004912040010008220180810000000000000008808000000000020400200800000200000000020000000000000000000000000000000000000000000000000000000000000000400800000000002000000000000000050000200000080000000000000000000020000800000000000008000000000000000000000000000000000000000000400000004000000000000000000000000000000000000004004000000000000000000000108008020000000000000000000000021");
      ByteFingerprint fp2 = new ByteFingerprint("fffffaff2b821f9880214a0088038100100100c004216008104000000200100000000000040010100400081e08007000070000000004068020010804004912040011000220080810000000000002008808040000000020400200800000200000000020008000004000000000000000000000000000000000000000000000000000000400800000000002000000000000000050000200000080000000000000000000020000800000000000008000000000000000000000000000000000000000000400000004000000000000000000000000000000000000004004000000000000000000000108008020000000000000000000000021");

      System.err.printf("n1=%d n2=%d sim=%f\n", fp1.getNBits(), fp2.getNBits(), fp1.tanimoto(fp2));
   }
}
