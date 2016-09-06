/*
   Copyright 2008-2014 Genentech Inc.

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
package com.genentech.chemistry.openEye;

import openeye.oechem.OEGraphMol;
import openeye.oechem.OEMolBase;

import com.aestel.chemistry.openEye.SimComparator;
import com.genentech.oechem.tools.OETools;

public abstract class AbstractAAPathComparatorTest
{  private String[] smis = new String[]
   {  "*", "C", "N", "CCO", "CC(=O)N", "c1ccccc1", "c1ncncc1", "c1[nH]ccc1",
      "c1ncncc1CC(=O)N", "c1ccccc1c1ncncc1"
   };

   private SimComparator<OEMolBase>[] molComps;
   private double[][] expectedSims;

   @SuppressWarnings("unchecked")
   public void setUp(double[][] expectedSims)
   {  this.expectedSims = expectedSims;
      // symetricize matrix
      for(int i=0; i<expectedSims.length; i++)
         for(int j=i+1; j<expectedSims.length; j++)
            expectedSims[j][i] = expectedSims[i][j];

      AAPathComparatorFact cFact = getComparatorFact();
      molComps = new SimComparator[smis.length];
      OEGraphMol mol = new OEGraphMol();
      for( int i=0; i<smis.length; i++)
      {  mol.Clear();
         OETools.smiToMol(mol, smis[i]);
         molComps[i] = cFact.createComparator(mol);
      }
      cFact.close();
   }


   abstract protected AAPathComparatorFact getComparatorFact();

   public void testComparator()
   {  // start with slightly more interesting case:
      int i1=1,j1=3;
      double sim1 = molComps[i1].similarity(molComps[j1]);
      assert Math.abs(sim1 - expectedSims[i1][j1]) < 0.001 :
         String.format("%2d,%2d\t%.4f != %.4f comparing %s\t%s\n",
                       i1,   j1,
                       sim1, expectedSims[i1][j1],
                       smis[i1], smis[j1]);

      for(int i=0; i< molComps.length; i++)
      {  for(int j=0; j<molComps.length; j++)
         {  double sim = molComps[i].similarity(molComps[j]);
            assert Math.abs(sim - expectedSims[i][j]) < 0.001 :
               String.format("%2d,%2d\t%.4f != %.4f comparing %s\t%s\n",
                             i,   j,
                             sim, expectedSims[i][j],
                             smis[i], smis[j]);
         }
      }
   }

   public void close()
   {  for( int i=0; i<molComps.length; i++)
         molComps[i].close();
   }
}
