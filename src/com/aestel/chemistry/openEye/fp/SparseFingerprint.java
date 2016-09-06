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

import java.util.*;

import com.aestel.utility.IntArrayList;

/**
 * Store fingerprint info by keeping the position of bits set to one.
 *
 * This is more efficient for spares fingerprints while {@link ByteFingerprint}
 * is more efficient for dense fingerprints.
 */
public class SparseFingerprint implements Fingerprint
{  /* sorted list of indexes of bits which are set. */
   private final int[] bits;

   public SparseFingerprint(int[] inBits)
   {  Arrays.sort(inBits);
      int last = Integer.MIN_VALUE;
      int countUnique = 0;
      for(int bit : inBits)
      {  if(bit != last)
         {  countUnique++;
            last = bit;
         }
      }

      last = Integer.MIN_VALUE;
      this.bits = new int[countUnique];
      countUnique = 0;
      for(int bit : inBits)
      {  if(bit != last)
         {  bits[countUnique++] = bit;
            last = bit;
         }
      }
   }

   /**
    * Conpute tanimoto by counting bits in common.
    */
   @Override
   public double tanimoto(Fingerprint other)
   {  if(other instanceof SparseFingerprint)
      {  int[] otherBits = ((SparseFingerprint)other).bits;
         int nBits1 = bits.length;
         int nBits2 = otherBits.length;

         int pos1 = 0;
         int pos2 = 0;
         int nCommon = 0;
         while(pos1 < bits.length && pos2 < otherBits.length)
         {  if(bits[pos1] == otherBits[pos2])
            {  nCommon++;
               pos1++;
               pos2++;
               continue;
            }
            if(bits[pos1] > otherBits[pos2])
               pos2++;
            else
               pos1++;
         }
         return ((double)nCommon) / (double) (nBits1 + nBits2 - nCommon);
      }

      throw new Error("not implemented yet");
   }

  @Override
public String getBinString()
   {  if(bits.length == 0) return "0";

      StringBuilder sb = new StringBuilder(bits[bits.length-1]+1);
      sb.setLength(bits[bits.length-1]+1);

      int bitIdx = 0;
      for(int pos=0; pos<sb.length(); pos++)
      {  if(pos < bits[bitIdx])
         {  sb.setCharAt(pos, '0');
         }else
         {  sb.setCharAt(pos, '1');
            bitIdx++;
         }
      }
      return sb.toString();
   }

   @Override
   public int[] getBits()
   {  return bits.clone();
   }

   /**
    * bit0 is on the far left.
    */
   @Override
   public String getHexString()
   {  if(bits.length == 0) return "0";

      StringBuilder sb = new StringBuilder((bits[bits.length-1]/4)+2);
      sb.setLength((bits[bits.length-1]/4)+1);
      int bitIdx = 0;
      int minNibbleBit = 0;    // smallest bit which fits into current nibble
      int bit = bits[bitIdx];
      for(int pos=0; pos<sb.length(); pos++)
      {  if(bit - minNibbleBit >= 4 )
         {  sb.setCharAt(pos, '0');
         }else
         {  int nibble = 0;
            while(bit - minNibbleBit < 4)
            {  nibble |= 1 << (minNibbleBit + 3 - bit);
               if(++bitIdx >= bits.length) break;
               bit = bits[bitIdx];
            }
            char nibbleChar = (char)(nibble < 10 ? '0' + nibble : 'a'-10+nibble);
            sb.setCharAt(pos, nibbleChar);
         }
         minNibbleBit += 4;
      }
      if(sb.length() % 2 == 1) sb.append('0');  // make even
      return sb.toString();
   }

   @Override
   public int getNBits()
   {  return bits.length;
   }

   /**
    * @return a folded with with newSize being the highest possible bit number.
    */
   @Override
   public Fingerprint fold(int newSize)
   {  IntArrayList newBits = new IntArrayList(bits.length);
      int firstShift = newSize  * 8 / 10;
      int secondShift = newSize * 64 / 100;
      int minFoldArea = newSize * 52 / 100;
      int sizeFoldArea= newSize - minFoldArea;
      for(int bit : bits)
      {  // do not fold onto first 20% of fingerprint because we assume these to
         // be the most frequent bits and therefore we assume these to be close
         // to saturated ( firstShift is size * 20%)
         if( bit >= newSize ) bit = bit - firstShift;
         // fold only once onto second 20% of newSize so that these do not get
         // too saturated
         if( bit >= newSize ) bit = bit - secondShift;
         // fold into area from minFoldArea to newSize
         if( bit >= newSize ) bit = minFoldArea + bit % sizeFoldArea;
         newBits.add(bit);
      }
      return new SparseFingerprint(newBits.toArray());
   }
}
