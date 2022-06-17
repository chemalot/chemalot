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

import openeye.oechem.OEExprOpts;
import openeye.oechem.OEGraphMol;
import openeye.oechem.OEMCSFunc;
import openeye.oechem.OEMCSSearch;
import openeye.oechem.OEMatchBase;
import openeye.oechem.OEMatchBaseIter;
import openeye.oechem.OEMolBase;
import openeye.oechem.oechem;

import com.aestel.chemistry.openEye.SimComparator;

/**
 * Returns the ratio of matched atoms to query atoms as similarity.
 * 
 * @author albertgo
 *
 */
public final class MCSSQueryRatioComparator implements SimComparator<OEMolBase>
{  private final OEMCSSearch mcss;
   private final OEMolBase pattern;
   private final int nAtomQuery;

   MCSSQueryRatioComparator(OEMolBase pattern, int type, int atExpr, int bdExpr, OEMCSFunc mcsFunc, int minAt)
   {  mcss = new OEMCSSearch(pattern, atExpr, bdExpr, type);
      this.pattern = new OEGraphMol(pattern);

      mcss.SetMCSFunc(mcsFunc);
      mcss.SetMinAtoms(minAt);
      nAtomQuery = pattern.NumAtoms();
   }

   @Override
   public double similarity(OEMolBase mol)
   {  int nAtomMatch = 0;

      OEMatchBaseIter mIt = mcss.Match(mol, true);
      if( mIt.hasNext() )
      {  OEMatchBase m = mIt.next();
         nAtomMatch = m.NumAtoms();
         m.delete();
      }

      mIt.delete();
      double sim = nAtomMatch / (double)nAtomQuery;

      return sim;
   }


   @Override
   public double similarity(SimComparator<OEMolBase> other)
   {  return similarity(((MCSSQueryRatioComparator) other).pattern);
   }

   /**
    * Default implementation to be improved.
    */
   @Override
   public double similarity(SimComparator<OEMolBase> other, double minSim)
   { return similarity(other);
   }

   @Override
   public void close()
   {  mcss.delete();
      pattern.delete();
   }


   public static void main(String argv[])
   {  String smiPat;
      String smiTar;

      smiPat = "c1cccn1C(=O)N1CCC1";
      smiTar = "c1cc[nH]c1";       // 0.194

      OEGraphMol mol1 = new OEGraphMol();
      oechem.OEParseSmiles(mol1, smiPat);
      System.err.println(mol1.IsValid() + " " + mol1.NumAtoms());

      OEGraphMol mol2 = new OEGraphMol();
      oechem.OEParseSmiles(mol2, smiTar);
      System.err.println(mol2.IsValid() + " " + mol2.NumAtoms());

      int atExpr = OEExprOpts.DefaultAtoms;
      int bdExpr = OEExprOpts.DefaultBonds;
      MCSSComparatorFact mcssCFact = new MCSSComparatorFact(MCSSCompareType.QUERYRatio, atExpr, bdExpr);
      SimComparator<OEMolBase> simer = mcssCFact.createComparator(mol1);
      try{  Thread.sleep(3000);
      } catch (InterruptedException e)
      {  e.printStackTrace();
      }
      long start = System.currentTimeMillis();
      double sim = simer.similarity(mol2);
      System.out.printf("%f\n", sim);
      for(int i=0; i<200; i++)
      {  // check if similarity is atom order independent
         oechem.OEScrambleMolecule(mol2);
         double sim2 = simer.similarity(mol2);
         if( sim2 != sim ) System.out.printf("Problem: sim is atom order dependent: %f\n", sim2);
      }
      System.err.printf("Done in %.3f\n", (System.currentTimeMillis()-start)/1000D);
      simer.close();
   }
}
