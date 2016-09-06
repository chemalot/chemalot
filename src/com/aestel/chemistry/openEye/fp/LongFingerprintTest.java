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

import org.testng.annotations.Test;

public class LongFingerprintTest

{
   @Test()
   public void testBin()
   {  assertBin(new long[] {0},                   "0000000000000000000000000000000000000000000000000000000000000000");
      assertBin(new long[] {0x8000000000000000L}, "1000000000000000000000000000000000000000000000000000000000000000");
      assertBin(new long[] {0xC000000000000000L}, "1100000000000000000000000000000000000000000000000000000000000000");
      assertBin(new long[] {0xC100000000000000L}, "1100000100000000000000000000000000000000000000000000000000000000");
      assertBin(new long[] {0x5000000000000000L}, "0101000000000000000000000000000000000000000000000000000000000000");
      assertBin(new long[] {0xffFFffFFffFFffFFL}, "1111111111111111111111111111111111111111111111111111111111111111");


      assertBin(new long[] {0,0x80},
          //|1234567|1234567|1234567|1234567|1234567|1234567|1234567|1234567|1234567|1234567|1234567|1234567|1234567|1234567|1234567|1234567
           "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000");
   }

   void assertBin(long[] bits, String res)
   {  LongFingerprint fp = new LongFingerprint(bits);
      assert fp.getBinString().equals(res)
         : "expected "+res+" got " + fp.getBinString();
   }

   @Test()
   public void testHex()
   {  assertHex(new long[] {0},                    "0000000000000000");
      assertHex(new long[] {0x8000000000000000L},  "8000000000000000");
      assertHex(new long[] {0xc000000000000000L},  "c000000000000000");
      assertHex(new long[] {0x9000000000000000L},  "9000000000000000");
      assertHex(new long[] {0x5000000000000000L},  "5000000000000000");
      assertHex(new long[] {0xd000000000000000L},  "d000000000000000");
      assertHex(new long[] {0x3000000000000000L},  "3000000000000000");
      assertHex(new long[] {0xb000000000000000L},  "b000000000000000");
      assertHex(new long[] {0x7000000000000000L},  "7000000000000000");
      assertHex(new long[] {0xf000000000000000L},  "f000000000000000");
      assertHex(new long[] {0x5200000000000000L},  "5200000000000000");
      assertHex(new long[] {0xd200000000000000L},  "d200000000000000");
      assertHex(new long[] {0x3200000000000000L},  "3200000000000000");
      assertHex(new long[] {0xb200000000000000L},  "b200000000000000");
      assertHex(new long[] {0x7200000000000000L},  "7200000000000000");
      assertHex(new long[] {0xf200000000000000L},  "f200000000000000");
      assertHex(new long[] {0x0800000000000000L},  "0800000000000000");
      assertHex(new long[] {0x4800000000000000L},  "4800000000000000");
      assertHex(new long[] {0x0c00000000000000L},  "0c00000000000000");
      assertHex(new long[] {0x0f00000000000000L},  "0f00000000000000");
      assertHex(new long[] {0,1},                  "00000000000000000000000000000001");
      assertHex(new long[] {0,0x8000000000000000L},"00000000000000008000000000000000");
      assertHex(new long[] {0x8000000000000000L,1},"80000000000000000000000000000001");
   }

   void assertHex(long[] bits, String res)
   {  LongFingerprint fp = new LongFingerprint(bits);
      assert fp.getHexString().equals(res)
         : "expected "+res+" got " + fp.getHexString();
   }

   @Test()
   public void testHex2()
   {  assertHex2(new LongFingerprint("0"),                 "0000000000000000");
      assertHex2(new LongFingerprint("7"),                 "7000000000000000");
      assertHex2(new LongFingerprint("700"),               "7000000000000000");
      assertHex2(new LongFingerprint("8000000000000000"),  "8000000000000000");
      assertHex2(new LongFingerprint("c000000000000000"),  "c000000000000000");
      assertHex2(new LongFingerprint("9000000000000000"),  "9000000000000000");
      assertHex2(new LongFingerprint("5000000000000000"),  "5000000000000000");
      assertHex2(new LongFingerprint("d000000000000000"),  "d000000000000000");
      assertHex2(new LongFingerprint("3000000000000000"),  "3000000000000000");
      assertHex2(new LongFingerprint("b000000000000000"),  "b000000000000000");
      assertHex2(new LongFingerprint("7000000000000000"),  "7000000000000000");
      assertHex2(new LongFingerprint("f000000000000000"),  "f000000000000000");
      assertHex2(new LongFingerprint("5200000000000000"),  "5200000000000000");
      assertHex2(new LongFingerprint("d220000000000000"),  "d220000000000000");
      assertHex2(new LongFingerprint("3220000000000000"),  "3220000000000000");
      assertHex2(new LongFingerprint("b220000000000000"),  "b220000000000000");
      assertHex2(new LongFingerprint("7220000000000000"),  "7220000000000000");
      assertHex2(new LongFingerprint("f220000000000000"),  "f220000000000000");
      assertHex2(new LongFingerprint("0880000000000000"),  "0880000000000000");
      assertHex2(new LongFingerprint("4880000000000000"),  "4880000000000000");
      assertHex2(new LongFingerprint("0cc0000000000000"),  "0cc0000000000000");
      assertHex2(new LongFingerprint("0ff0000000000000"),  "0ff0000000000000");
      assertHex2(new LongFingerprint("00000000000000001"), "00000000000000001000000000000000");
      assertHex2(new LongFingerprint("000000000000000080"),"00000000000000008000000000000000");
   }

   void assertHex2(LongFingerprint fp, String res)
   {  assert fp.getHexString().equals(res)
         : "expected "+res+" got " + fp.getHexString();
   }

   @Test()
   public void testTanimoto()
   {  assert new LongFingerprint("7F8")
         .tanimoto(new LongFingerprint("7F8"))
         == 1D;

      assert new LongFingerprint("FFA")        //11111111 11
         .tanimoto(new LongFingerprint("7FA")) //01111111 11
         == 9D/10D;

      assert new LongFingerprint("FFA")         //1111 1111 1010
         .mtanimoto(new LongFingerprint("7FA")) //0111 1111 1010
         == 9D/11D;

      assert new LongFingerprint("FFA1")       //1111 1111 1010 0001
         .tanimoto(new LongFingerprint("7FA")) //0111 1111 1010 0000
         == 9D/(11+9-9);

      // ed80000040000000800000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001
      // edc000004000000080
      assert new LongFingerprint(
            new SparseFingerprint(new int[] {0,1,2,33,4,5,64,7,8,999}).getHexString())
         .tanimoto(new LongFingerprint(
               new SparseFingerprint(new int[] { 0,1,2,33,4,5,64,7,8,9}).getHexString()))
         == 9D/11D;
   }

//   @Test
//   public void testFold()
//   {  assertFold(new int [] {0},    "80");
//      assertFold(new int [] {102},  "00000000000000000000000002");
//      assertFold(new int [] {103},  "00000000000000000000000001");
//      assertFold(new int [] {104},  "0000000000000000000000000080");
//      assertFold(new int [] {511},  "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001");
//      assertFold(new int [] {512},  "00000000000000000000000001");
//      assertFold(new int [] {513},  "0000000000000000000000000080");
//      assertFold(new int [] {920},  "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001");
//      assertFold(new int [] {921},  "000000000000000000000000000000000000000000000040");
//      assertFold(new int [] {922},  "000000000000000000000000000000000000000000000020");
//      assertFold(new int [] {1247}, "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001");
//      assertFold(new int [] {1248}, "000000000000000000000000000000000000000000000000000000000000000000000002");
//      assertFold(new int [] {1249}, "000000000000000000000000000000000000000000000000000000000000000000000001");
//      assertFold(new int [] {1472}, "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000002");
//      assertFold(new int [] {1473}, "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001");
//      assertFold(new int [] {1474}, "00000000000000000000000000000000000000000000000000000000000000000020");
//   }

//   private void assertFold(int[] bits, String expectedHex)
//   {  SparseFingerprint fp = new SparseFingerprint(bits);
//      String res = fp.fold(512).getHexString();
//
//      assert expectedHex.equals(res)
//         : "Expected: " + expectedHex + " got: " + res;
//   }


   @Test()
   public void testCount()
   {  assert new LongFingerprint("FFA").getNBits()== 10;
      assert new LongFingerprint("FFA").getNBits()== 10;
      assert new LongFingerprint("00008001080").getNBits()== 3;
   }
}
