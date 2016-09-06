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
import openeye.oechem.OEMCSMaxAtoms;
import openeye.oechem.OEMCSSearch;
import openeye.oechem.OEMCSType;
import openeye.oechem.OEMatchBase;
import openeye.oechem.OEMatchBaseIter;
import openeye.oechem.oechem;

public final class MCSSSimilarity
{  private final OEMCSSearch mcss;
   private final int nAtomQuery;
   private final int nBondQuery;

   public static MCSSSimilarity factory(OEGraphMol mol, String typeName)
   {  if( "DEFAULT".equalsIgnoreCase(typeName) )
      {  int type = OEMCSType.Approximate;
         int minAt = 2;
         int atExpr = OEExprOpts.DefaultAtoms;
         int bdExpr = OEExprOpts.DefaultBonds;
         OEMCSMaxAtoms mcsFunc = new OEMCSMaxAtoms();


         return new MCSSSimilarity(mol, type, atExpr, bdExpr, mcsFunc, minAt);
      }else
      {  throw new Error("Unknown type: " +typeName);
      }
   }

   private MCSSSimilarity(OEGraphMol pattern, int type, int atExpr, int bdExpr, OEMCSFunc mcsFunc, int minAt)
   {  mcss = new OEMCSSearch(pattern, atExpr, bdExpr, type);
      mcss.SetMCSFunc(mcsFunc);
      mcss.SetMinAtoms(minAt);
//System.err.println(mcss.IsValid());
      nAtomQuery = pattern.NumAtoms();
      nBondQuery = pattern.NumBonds();
   }

   public double similarity(OEGraphMol mol)
   {  int nAtomMatch = 0;
      int nBondMatch = 0;

      OEMatchBaseIter mIt = mcss.Match(mol, true);
      if( mIt.hasNext() )
      {  OEMatchBase m = mIt.next();
         nAtomMatch = m.NumAtoms();
         nBondMatch = m.NumBonds();
         m.delete();
      }
      mIt.delete();
      double sim = nAtomMatch / (double)(nAtomQuery + mol.NumAtoms() - nAtomMatch);
      sim = sim + nBondMatch / (double)(nBondQuery + mol.NumBonds() - nBondMatch);
      sim /= 2;

      return sim;
   }

   public void close()
   {  mcss.delete();
   }


   public static void main(String argv[])
   {  String smiPat;
      String smiTar;

      // aa paper 1
      smiPat = "c1ccn[nH]1";
      smiTar = "Cc1ccno1";

      // aa paper 2
      smiPat = "c1ccnn1C(=O)N1CCC1";
      smiTar = "c1cc[nH]c1C(=O)N1CCC1"; // 0.321
      smiTar = "c1cc[nH]c1";       // 0.057

      smiPat = "c1cccn1C(=O)N1CCC1";
      smiTar = "c1cc[nH]c1";       // 0.194

      OEGraphMol pattern = new OEGraphMol();
      oechem.OEParseSmiles(pattern, smiPat);
      System.err.println(pattern.IsValid() + " " + pattern.NumAtoms());

      OEGraphMol target = new OEGraphMol();
      oechem.OEParseSmiles(target, smiTar);
      System.err.println(target.IsValid() + " " + target.NumAtoms());

      MCSSSimilarity simer = MCSSSimilarity.factory(pattern, "DEFAULT");
      System.out.printf("%f\n", simer.similarity(target));
      simer.close();
   }


}
