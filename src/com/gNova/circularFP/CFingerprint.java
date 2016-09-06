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

package com.gNova.circularFP;

import java.util.*;
import java.util.Map.Entry;

/**
 *
 * @author tjodonnell
 *
 */

public class CFingerprint implements Fingerprint
{  enum CFPCountType { NOCount, LOGCount, LINCount }
   private Set<Integer> fp;
   private int length;

   public static CFingerprint createCFingerprint(Map<Integer,Integer>fp, int length, CFPCountType countType)
   {  switch(countType)
      {  case NOCount:
            return new CFingerprint(fp.keySet(), length);
         case LINCount:
            return createLinCountFP(fp, length);
         case LOGCount:
            return createLogCountFP(fp, length);
         default:
            assert false;
      }
      return null;
   }

   private CFingerprint(Set<Integer>fp, int length)
   {  this.fp = new HashSet<Integer>(fp);
      this.length = length;
   }

   /**
    * Create fingerprint that also includes bits for replicates of bits in log2 steps
    * @param fp Map with atomIds and counts
    * @param length
    */
   private static CFingerprint createLogCountFP(Map<Integer,Integer>fp, int length)
   {  HashSet<Integer> countfp = new HashSet<Integer>(fp.size()*2);
      int mulitplier = Prime.getPrimeLargerEqThan(length*2+7);

      for(Entry<Integer, Integer> e : fp.entrySet())
      {  int count = e.getValue();
         int pos = e.getKey();
         while(count > 0)
         {  countfp.add(pos);
            pos = pos * mulitplier;    // compute new hash position
            count /= 2;                // set additional bits doubling bit counts
         }
      }
      return new CFingerprint(countfp,length);
   }

   /**
    * Create fingerprint that also includes bits for replicates of bits in linear steps
    * @param fp Map with atomIds and counts
    * @param length
    */
   private static CFingerprint createLinCountFP(Map<Integer,Integer>fp, int length)
   {  HashSet<Integer> countfp = new HashSet<Integer>(fp.size()*2);
      int mulitplier = Prime.getPrimeLargerEqThan(length*2+7);

      for(Entry<Integer, Integer> e : fp.entrySet())
      {  int count = e.getValue();
         int pos = e.getKey();
         while(count > 0)
         {  countfp.add(pos);
            pos = pos * mulitplier;    // compute new hash position
            count--;                   // set additional bits doubling bit counts
         }
      }
      return new CFingerprint(countfp,length);
   }

   public List<Integer> getBitNums()
   {
      List<Integer> nums = new ArrayList<Integer>();
      for (Integer bit : this.fp)
      {
         int b = (bit&Integer.MAX_VALUE) % this.length;
         nums.add(b);
      }

      // tests
//      nums.clear();
//      nums.add(0);
//      nums.add(31);
//      nums.add(32);
//      nums.add(this.length-1);

      return nums;
   }


   @Override
   public int[] getBits()
   {
      List<Integer> nums = this.getBitNums();
      Collections.sort(nums);
      int[] x = new int[nums.size()];
      int i=0;
      for (Integer bit : nums)
      {
         x[i++] = bit;
      }
      return x;
   }

   @Override
   public BitSet getBitSet()
   {
      BitSet bset = new BitSet(this.length);
      for (Integer bit : this.getBitNums())
      {
         bset.set(bit);
      }
      return bset;
   }

   @Override
   public String getBitString()
   {
      BitSet nums = this.getBitSet();
      //Collections.sort(nums);
      String x = nums.toString().replace(" ","");
      return x.substring(1, x.length()-1);
   }

   @Override
   public String getHexString()
   {
      // bit 0 on far left
      int nelem = (int)Math.ceil(this.length/32.);
      int[] a = new int[nelem];
      for (Integer bit : this.getBitNums())
      {
         a[bit/32] |=  1 << (31 - bit % 32);
      }
      String hex = "";
      for (int i=0; i<nelem; ++i)
      {
         hex += String.format("%08x", a[i]);
      }
      return hex;
   }

   @Override
   public String getAtomIDString()
   {
      String x = this.fp.toString().replace(" ","");
      return x.substring(1, x.length()-1);
   }

   @Override
   public double Tanimoto(Fingerprint other)
   {
      return this.Tanimoto((CFingerprint)other);
   }

   public double Tanimoto(CFingerprint other)
   {
      if (this.length != other.length) return -1.0;
      BitSet A;
      BitSet B;

      A = this.getBitSet();
      B = other.getBitSet();
      A.and(B); // modifies A
      double c = A.cardinality();

      A = this.getBitSet();
      B = other.getBitSet();
      B.andNot(A); // modifies B
      double b = B.cardinality();

      A = this.getBitSet();
      B = other.getBitSet();
      A.andNot(B); // modifies A
      double a = A.cardinality();

      return c / (a + b + c);
   }

   @Override
   public int getNBits()
   {
      return this.fp.size();
   }

}
