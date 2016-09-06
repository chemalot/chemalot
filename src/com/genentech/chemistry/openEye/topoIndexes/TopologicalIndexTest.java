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
package com.genentech.chemistry.openEye.topoIndexes;

import java.io.IOException;
import java.util.concurrent.TimeoutException;

import openeye.oechem.OEGraphMol;
import openeye.oechem.OEMolBase;
import openeye.oechem.oechem;

import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import com.genentech.oechem.tools.OETools;

public class TopologicalIndexTest
{
   private TopologicalIndexer tp;
   private OEGraphMol mol;

   @BeforeClass
   public void setUp() throws IOException
   {  tp = new TopologicalIndexer();
      mol = new OEGraphMol();
   }

   @Test
   public void testNPentan() throws TimeoutException
   {  testSmi("[Ag]BrC[D]F", 2.1906, 1.9283, 3.7619, 2.1906, 1.9283, 3.7619, 20, 14);

      testSmi("CCCCC", 2.1906, 2.1906, 2.1906, 2.1906, 2.1906, 2.1906, 20, 14);

   }

   @Test
   public void testIButan() throws TimeoutException
   {  testSmi("[Ag]BrC([D])F", 2.539539, 2.2677, 4.1403,  2.539539, 2.2677, 4.1403,  18, 16);
   }

   @Test
   public void testCycPentan() throws TimeoutException
   {  testSmi("[Ag]1BrC[D]F1", 2.083333, 1.78116, 4.5139,  0.7813, 0.6679, 1.6927,  15, 20);
   }

   @Test
   public void testMetCycPentan() throws TimeoutException
   {  testSmi("[Ag]1BrC[D]F1[Ge]", 2.184106, 1.8837, 4.3537,  0.9051, 0.7809, 1.7606,  26, 26);
   }

   @Test
   public void test3MetCycPenten() throws TimeoutException
   {  testSmi("[Ag]1=BrC[D]F1[Ge]", 2.377489, 2.0379, 4.8359,  0.996, 0.8602, 1.9339,  26, 26);
   }

   @Test
   public void test1MetCycPenten() throws TimeoutException
   {  testSmi("[Ag]1BrC[D]=F1[Ge]", 2.4496, 2.1265, 4.7972,  0.9964, 0.8588, 1.9451,  26, 26);
   }

   @Test
   public void test1MetCycPentenO() throws TimeoutException
   {  testSmi("C1=CC(C)CO1", 2.3775, 2.4069, 2.2872, 1.0095,  1.0223, 0.9697,  26, 26);
   }

   @Test
   public void testNonHetero() throws TimeoutException
   {  testSmi("CCCCCCCCCCCCC(C)C1=CCC=C1CC",  1.9356, 1.9356, 1.9356,  1.6812, 1.6807, 1.6818,  1303,  90);
      testSmi("CCCCCCCCCCCCC(C)C1=CNC=C1C#N", 1.9820, 1.9870, 1.9646,  1.7512, 1.7548, 1.7376,  1303,  90);
   }

   @Test
   public void testPyrol() throws TimeoutException
   {  testSmi("c1c[nH]cc1", 3.125, 3.1491, 3.0414, 1.1719, 1.1806, 1.1409, 15, 20);
      testSmi("C1=CNC=C1",  3.125, 3.1491, 3.0414, 1.1719, 1.1806, 1.1409, 15, 20);
   }

   @Test
   public void testHetero() throws TimeoutException
   {  testSmi("FBCNO[Si]Cl[As][Se]Br[Te]IC(=O)C1=CNC=C1C#N", 1.9958, 1.8266, 1.8766,  1.7623, 1.6082, 1.6548,  1303, 90);
   }

   @Test
   public void testoMetToluol() throws TimeoutException
   {  testSmi("Cc1cc(C)ccc1", 3.0776, 3.0776, 3.0789,  1.4467, 1.4462, 1.4472,  61, 36);
   }

   @Test
   public void testFF() throws TimeoutException
   {  testSmi("FF",       1, 1.1271, 0.6914,  1, 1.1271, 0.6914,  1, 2);
   }

   @Test
   public void testBB() throws TimeoutException
   {   testSmi("BB",       1, 0.9634, 1.1755,  1, 0.9634, 1.1755,  1, 2);
   }

   @Test
   public void testCC() throws TimeoutException
   {   testSmi("CC",       1, 1, 1,   1, 1, 1,      1, 2);
   }

   @Test
   public void testNN() throws TimeoutException
   {   testSmi("NN",       1, 1.0389, 0.8701, 1, 1.0389, 0.8701,  1, 2);
   }

   @Test
   public void testOO() throws TimeoutException
   {   testSmi("OO",       1, 1.0812, 0.7708,  1, 1.0812, 0.7708,  1, 2);
   }

   @Test
   public void testSiSi() throws TimeoutException
   {   testSmi("[Si][Si]", 1, 0.8862, 1.0670,  1, 0.8862, 1.0670,  1, 2);
   }

   @Test
   public void testClCl() throws TimeoutException
   {   testSmi("ClCl",     1, 0.9859, 0.7226,  1, 0.9859, 0.7226,  1, 2);
   }

   @Test
   public void testAsAs() throws TimeoutException
   {   testSmi("[As][As]", 1, 0.7254, 1.0574,  1, 0.7254, 1.0574,  1, 2);
   }

   @Test
   public void testSeSe() throws TimeoutException
   {   testSmi("[Se][Se]", 1, 0.7458, 0.9136,  1, 0.7458, 0.9136,  1, 2);
   }

   @Test
   public void testBrBr() throws TimeoutException
   {   testSmi("BrBr",     1, 0.7673, 0.8042,  1, 0.7673, 0.8042,  1, 2);
   }

   @Test
   public void testTeTe() throws TimeoutException
   {   testSmi("[Te][Te]", 1, 0.6139, 1.0480,  1, 0.6139, 1.0480,  1, 2);
   }

   @Test
   public void testII() throws TimeoutException
   {   testSmi("II",       1, 0.6285, 0.9065,  1, 0.6285, 0.9065,  1, 2);
}

   @Test
   public void testCdO() throws TimeoutException
   {   testSmi("C=O",      2, 2.0796, 1.7559,  2, 2.0796, 1.7559,  1, 2);
   }

   @Test
   public void testCtN() throws TimeoutException
   {   testSmi("C#N",      3, 3.0577, 2.7993,  3, 3.0577, 2.7993,  1, 2);
   }

   @Test
   static void testCache() throws TimeoutException
   {  OEMolBase mol = new OEGraphMol();//.Clear();
      oechem.OEParseSmiles(mol,"C#N");
      TopologicalIndexer tp2 = new TopologicalIndexer();
      tp2.computeIndexes(mol);
      Assert.assertEquals(tp2.getBalabanJIndex(),  3, 0.0001D, "J "  + oechem.OECreateSmiString(mol));
      tp2.close();

      tp2 = new TopologicalIndexer();
      tp2.computeIndexes(mol);
      Assert.assertEquals(tp2.getBalabanJXIndex(),  2.7993, 0.001D, "JX "  + oechem.OECreateSmiString(mol));
      tp2.close();

      tp2 = new TopologicalIndexer();
      tp2.computeIndexes(mol);
      Assert.assertEquals(tp2.getBalabanJYIndex(),  3.0577, 0.001D, "JY "  + oechem.OECreateSmiString(mol));
      tp2.close();

      tp2 = new TopologicalIndexer();
      tp2.computeIndexes(mol);
      Assert.assertEquals(tp2.getBalabanJStarIndex(),  3, 0.0001D, "J* "  + oechem.OECreateSmiString(mol));
      tp2.close();

      tp2 = new TopologicalIndexer();
      tp2.computeIndexes(mol);
      Assert.assertEquals(tp2.getBalabanJXStarIndex(),  2.7993, 0.001D, "JX* "  + oechem.OECreateSmiString(mol));
      tp2.close();

      tp2 = new TopologicalIndexer();
      tp2.computeIndexes(mol);
      Assert.assertEquals(tp2.getBalabanJYStarIndex(),  3.0577, 0.001D, "JY* "  + oechem.OECreateSmiString(mol));
      tp2.close();

      tp2 = new TopologicalIndexer();
      tp2.computeIndexes(mol);
      Assert.assertEquals(tp2.getWienerIndex(),  1, "Wiener "  + oechem.OECreateSmiString(mol));
      tp2.close();

      tp2 = new TopologicalIndexer();
      tp2.computeIndexes(mol);
      Assert.assertEquals(tp2.getZagrebIndex(),  2, "ZagrebM1 "  + oechem.OECreateSmiString(mol));
      tp2.close();

      mol.delete();
   }

   private void testSmi(String smi, double expectedJ, double expectedJY, double expectedJX,
            double expectedJSt, double expectedJYSt, double expectedJXSt,
            int expectedWiener, int expectedZagreb) throws TimeoutException
   {  mol.Clear();
      oechem.OEParseSmiles(mol,smi);
      tp.computeIndexes(mol);
      System.err.printf("%-25s BJ=%7.4f BJY=%7.4f BJX=%7.4f BJ*=%7.4f BJY*=%7.4f BJX*=%7.4f w=%4d m1=%4d\n", OETools.molToCanSmi(mol, true),
               tp.getBalabanJIndex(),     tp.getBalabanJYIndex(),     tp.getBalabanJXIndex(),
               tp.getBalabanJStarIndex(), tp.getBalabanJYStarIndex(), tp.getBalabanJXStarIndex(),
               tp.getWienerIndex(), tp.getZagrebIndex());
      Assert.assertEquals(tp.getBalabanJIndex(),      expectedJ,    0.0001D, "J " + oechem.OECreateSmiString(mol));
      Assert.assertEquals(tp.getBalabanJYIndex(),     expectedJY,   0.001D, "JY " + oechem.OECreateSmiString(mol));
      Assert.assertEquals(tp.getBalabanJXIndex(),     expectedJX,   0.001D, "JX " + oechem.OECreateSmiString(mol));
      Assert.assertEquals(tp.getBalabanJStarIndex(),  expectedJSt,  0.0001D, "J " + oechem.OECreateSmiString(mol));
      Assert.assertEquals(tp.getBalabanJYStarIndex(), expectedJYSt, 0.001D, "JY " + oechem.OECreateSmiString(mol));
      Assert.assertEquals(tp.getBalabanJXStarIndex(), expectedJXSt, 0.001D, "JX " + oechem.OECreateSmiString(mol));
      Assert.assertEquals(tp.getWienerIndex(), expectedWiener, "Wiener " + oechem.OECreateSmiString(mol));
      Assert.assertEquals(tp.getZagrebIndex(), expectedZagreb, "ZagrebM1" + oechem.OECreateSmiString(mol));
   }

   @AfterClass
   public void close()
   {  mol.delete();
      mol = null;

      tp.close();
      tp = null;
   }
}
