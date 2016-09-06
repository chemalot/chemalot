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

public class ByteFingerprintTest

{
   @Test()
   public void testBin()
   {  assertBin(new byte[] {0},            "00000000");
      assertBin(new byte[] {(byte)0x80},   "10000000");
      assertBin(new byte[] {(byte)0xC0},   "11000000");
      assertBin(new byte[] {(byte)0xC1},   "11000001");
      assertBin(new byte[] {(byte)0x50},   "01010000");


      assertBin(new byte[] {0,0,0,0,0,0,0,0,(byte)0x80},
          //|1234567|1234567|1234567|1234567|1234567|1234567|1234567|1234567|1234567
           "000000000000000000000000000000000000000000000000000000000000000010000000");
   }

   void assertBin(byte[] bits, String res)
   {  ByteFingerprint fp = new ByteFingerprint(bits);
      assert fp.getBinString().equals(res)
         : "expected "+res+" got " + fp.getBinString();
   }

   @Test()
   public void testHex()
   {  assertHex(new byte[] {0},                         "00");
      assertHex(new byte[] {(byte)0x80},                "80");
      assertHex(new byte[] {(byte)0xc0},                "c0");
      assertHex(new byte[] {(byte)0x90},                "90");
      assertHex(new byte[] {(byte)0x50},                "50");
      assertHex(new byte[] {(byte)0xd0},                "d0");
      assertHex(new byte[] {(byte)0x30},                "30");
      assertHex(new byte[] {(byte)0xb0},                "b0");
      assertHex(new byte[] {(byte)0x70},                "70");
      assertHex(new byte[] {(byte)0xf0},                "f0");
      assertHex(new byte[] {(byte)0x52},                "52");
      assertHex(new byte[] {(byte)0xd2},                "d2");
      assertHex(new byte[] {(byte)0x32},                "32");
      assertHex(new byte[] {(byte)0xb2},                "b2");
      assertHex(new byte[] {(byte)0x72},                "72");
      assertHex(new byte[] {(byte)0xf2},                "f2");
      assertHex(new byte[] {(byte)0x08},                "08");
      assertHex(new byte[] {(byte)0x48},                "48");
      assertHex(new byte[] {(byte)0x0c},                "0c");
      assertHex(new byte[] {(byte)0x0f},                "0f");
      assertHex(new byte[] {0,0,0,0,0,0,0,1},           "0000000000000001");
      assertHex(new byte[] {0,0,0,0,0,0,0,0,(byte)0x80},"000000000000000080");
   }

   void assertHex(byte[] bits, String res)
   {  ByteFingerprint fp = new ByteFingerprint(bits);
      assert fp.getHexString().equals(res)
         : "expected "+res+" got " + fp.getHexString();
   }

   @Test()
   public void testHex2()
   {  assertHex2(new ByteFingerprint("0"),                  "00");
      assertHex2(new ByteFingerprint("7"),                  "70");
      assertHex2(new ByteFingerprint("70"),                 "70");
      assertHex2(new ByteFingerprint("700"),              "7000");
      assertHex2(new ByteFingerprint("80"),                 "80");
      assertHex2(new ByteFingerprint("c0"),                 "c0");
      assertHex2(new ByteFingerprint("90"),                 "90");
      assertHex2(new ByteFingerprint("50"),                 "50");
      assertHex2(new ByteFingerprint("d0"),                 "d0");
      assertHex2(new ByteFingerprint("30"),                 "30");
      assertHex2(new ByteFingerprint("b0"),                 "b0");
      assertHex2(new ByteFingerprint("70"),                 "70");
      assertHex2(new ByteFingerprint("f0"),                 "f0");
      assertHex2(new ByteFingerprint("52"),                 "52");
      assertHex2(new ByteFingerprint("d2"),                 "d2");
      assertHex2(new ByteFingerprint("32"),                 "32");
      assertHex2(new ByteFingerprint("b2"),                 "b2");
      assertHex2(new ByteFingerprint("72"),                 "72");
      assertHex2(new ByteFingerprint("f2"),                 "f2");
      assertHex2(new ByteFingerprint("08"),                 "08");
      assertHex2(new ByteFingerprint("48"),                 "48");
      assertHex2(new ByteFingerprint("0c"),                 "0c");
      assertHex2(new ByteFingerprint("0f"),                 "0f");
      assertHex2(new ByteFingerprint("0000000000000001"),   "0000000000000001");
      assertHex2(new ByteFingerprint("000000000000000080"), "000000000000000080");
   }

   void assertHex2(ByteFingerprint fp, String res)
   {  assert fp.getHexString().equals(res)
         : "expected "+res+" got " + fp.getHexString();
   }

   @Test()
   public void testTanimoto()
   {  assert new ByteFingerprint("7F8")
         .tanimoto(new ByteFingerprint("7F8"))
         == 1D;

      assert new ByteFingerprint("FFA")        //11111111 11
         .tanimoto(new ByteFingerprint("7FA")) //01111111 11
         == 9D/10D;

      assert new ByteFingerprint("FFA")         //1111 1111 1010 0000
         .mtanimoto(new ByteFingerprint("7FA")) //0111 1111 1010 0000
         == 9D/11D;

      assert new ByteFingerprint("FFA1")      // 1111 1111 1010 0001
         .tanimoto(new ByteFingerprint("7FA"))// 0111 1111 1010 0000
         == 9D/11D;

      assert new ByteFingerprint(
            new SparseFingerprint(new int[] {0,1,2,33,4,5,64,7,8,999}).getHexString())
         .tanimoto(new ByteFingerprint(
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
   {  assert new ByteFingerprint("FFA").getNBits()== 10;
      assert new ByteFingerprint("FFA").getNBits()== 10;
      assert new ByteFingerprint("00008001080").getNBits()== 3;
   }
}
