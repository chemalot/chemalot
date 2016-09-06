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

public class SparseFingerprintTest

{  
   @Test()
   public void testBin()
   {  assertBin(new int[] {},      "0");
      assertBin(new int[] {0},     "1");
      assertBin(new int[] {0,1},   "11");
      assertBin(new int[] {0,1,8}, "110000001");
      assertBin(new int[] {4,1},   "01001");
      assertBin(new int[] {65},    "000000000000000000000000000000000000000000000000000000000000000001");
   }
 
   void assertBin(int[] bits, String res)
   {  SparseFingerprint sp = new SparseFingerprint(bits);
      assert sp.getBinString().equals(res) 
         : "expected "+res+" got " + sp.getBinString();
   }
   
   @Test()
   public void testHex()
   {  assertHex(new int[] {},         "0");
      assertHex(new int[] {0},        "80");
      assertHex(new int[] {0,1},      "c0");
      assertHex(new int[] {0,3},      "90");
      assertHex(new int[] {1,3},      "50");
      assertHex(new int[] {0,1,3},    "d0");
      assertHex(new int[] {2,3},      "30");
      assertHex(new int[] {0,2,3},    "b0");
      assertHex(new int[] {1,2,3},    "70");
      assertHex(new int[] {0,1,2,3},  "f0");
      assertHex(new int[] {6,1,3},    "52");
      assertHex(new int[] {6,0,1,3},  "d2");
      assertHex(new int[] {6,2,3},    "32");
      assertHex(new int[] {6,0,2,3},  "b2");
      assertHex(new int[] {6,1,2,3},  "72");
      assertHex(new int[] {6,0,1,2,3},"f2");
      assertHex(new int[] {4},        "08");
      assertHex(new int[] {4,1},      "48");
      assertHex(new int[] {4,5},      "0c");
      assertHex(new int[] {4,5,6,7},  "0f");
      assertHex(new int[] {63},       "0000000000000001");
      assertHex(new int[] {64},       "000000000000000080");
   }

   void assertHex(int[] bits, String res)
   {  SparseFingerprint sp = new SparseFingerprint(bits);
      assert sp.getHexString().equals(res) 
         : "expected "+res+" got " + sp.getHexString();
   }

   @Test()
   public void testTanimoto()
   {  assert new SparseFingerprint(new int[] {0,1,2,3,4,5,6,7,8,9})
         .tanimoto(new SparseFingerprint(new int[] { 0,1,2,3,4,5,6,7,8,9}))
         == 1D;
   
      assert new SparseFingerprint(new int[] {0,1,2,3,4,5,6,7,8,9})
         .tanimoto(new SparseFingerprint(new int[] { 1,2,3,4,5,6,7,8,9}))
         == .9D;

      assert new SparseFingerprint(new int[] {0,1,2,33,4,5,64,7,8,999})
         .tanimoto(new SparseFingerprint(new int[] { 0,1,2,33,4,5,64,7,8,9}))
         == 9D/11D;
   }
   
   @Test
   public void testFold()
   {  assertFold(new int [] {0},    "80");
      assertFold(new int [] {102},  "00000000000000000000000002");
      assertFold(new int [] {103},  "00000000000000000000000001");
      assertFold(new int [] {104},  "0000000000000000000000000080");
      assertFold(new int [] {511},  "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001");
      assertFold(new int [] {512},  "00000000000000000000000001");
      assertFold(new int [] {513},  "0000000000000000000000000080");
      assertFold(new int [] {920},  "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001");
      assertFold(new int [] {921},  "000000000000000000000000000000000000000000000040");
      assertFold(new int [] {922},  "000000000000000000000000000000000000000000000020");
      assertFold(new int [] {1247}, "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001");
      assertFold(new int [] {1248}, "000000000000000000000000000000000000000000000000000000000000000000000002");
      assertFold(new int [] {1249}, "000000000000000000000000000000000000000000000000000000000000000000000001");
      assertFold(new int [] {1472}, "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000002");
      assertFold(new int [] {1473}, "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001");
      assertFold(new int [] {1474}, "00000000000000000000000000000000000000000000000000000000000000000020");
   }
   
   private void assertFold(int[] bits, String expectedHex)
   {  SparseFingerprint fp = new SparseFingerprint(bits);
      String res = fp.fold(512).getHexString();
      
      assert expectedHex.equals(res)
         : "Expected: " + expectedHex + " got: " + res;
   }
   
   
   @Test()
   public void testCount()
   {  assert new SparseFingerprint(new int[] {0,1,2,3,4,5,6,7,8,9}).getNBits() == 10;
      assert new SparseFingerprint(new int[] {0,0,1111111,0,0,0,1}).getNBits() == 3;
      assert new SparseFingerprint(new int[] {0,0,0,0,0,0,0,0,0,0}).getNBits() == 1;
   }
}
