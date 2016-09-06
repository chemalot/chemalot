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
import openeye.oechem.oechem;

import com.aestel.chemistry.openEye.SimComparator;
import com.aestel.math.HungarianAlgorithm;
import com.genentech.chemistry.openEye.AAPathComparatorFact.AAPathCompareType;

/**
 * Compute Atom Atom Path similarity.
 *
 * First a similarity matrix is computed giving a similarity of each atom in the
 * Comparators molecule with each atom of the database molecule. The similarity
 * is defined as the tanimoto coefficient of the paths the atoms share/not share.
 *
 * Then each atom of the smaller molecule is matched to one atom of the other
 * based on maximum atom pair similarity.
 *
 * This implementation use the hungarian algorithm to find the optimal assignment
 * of atoms and to compute the similarity of two atoms is changed as follows (same as comparator2):
 *    pathsInCommon / (2 * max(nPath(AtomMolA), nPath(AtomMolB)) - pathsInCommon)
 *
 *
 * The total similarity is computed as
 *    sum(Matched(atomPairSimilarity))/ (2 * max(nAtomsMolA, nAtomsMolB) - sum(Matched(atomPairSimilarity)))
 *
 * @author albertgo 2012
 *
 */
public class AAPathComparator4 extends AAPathComparator2
{  @SuppressWarnings("unused")
   private static final boolean DEBUG = false;

   AAPathComparator4(OEMolBase mol, int maxBonds, double minAtSim)
   {  super(mol, maxBonds, minAtSim);
   }


   @Override
   protected double similarity(int nAtoms1, String[][] atPath1, int nAtoms2,
            String[][] atPath2)
   { // make user atPath length is smaller than atPath2 length
      if (nAtoms1 > nAtoms2)
      {  String[][] d = atPath1;
         atPath1 = atPath2;
         atPath2 = d;

         int nd = nAtoms1;
         nAtoms1 = nAtoms2;
         nAtoms2 = nd;
//         System.err.println("Swapped");
      }

      double[][] atSims = getAtomSimilarities(atPath1, atPath2);

      // the oechem toolkit can have atom indexes which do not have atoms
      // ie. mol.GetMaxAtomIdx() returns more than mol.NumAtoms()
      // for the hungarian algorithm we need a matrix which :
      // - has one row per atoms in nAtoms1
      // - has one column per atoms in nAtoms2
      // - is filled up with high costs to make a nAtoms2 x nAtoms2 matrix

      // create square matrix with distances as hungarian algorithm needs
      // square matrix and it needs to minimize the cost (=distance)
      double[][] distMat;

      distMat = new double[nAtoms2][nAtoms2];
      // fill in non-square part of the distance matrix with costs 2 and up
      // this will make the atoms not existing in atPath1 be asigned last
      for (int i = nAtoms1; i < nAtoms2; i++)
         for (int j = 0; j < nAtoms2; j++)
            distMat[i][j] = 2D; // higher cost than similarity=0 (dist=1)

      // map from indexes 0-nAtomsX indexes to
      // 0-atPathX.length which may have non-existing atoms to
      int at1IdxMap[] = new int[nAtoms1];
      int at2IdxMap[] = new int[nAtoms2];

      // convert similarities to distance
      for (int i = 0, iAt = 0; i < atPath1.length; i++)
      {  if (atPath1[i] != null) // deal with empty atoms in oechem
         {  for (int j = 0, jAt = 0; j < atPath2.length; j++)
            {  if (atPath2[j] != null)
               {  distMat[iAt][jAt] = 1 - atSims[i][j];
                  at2IdxMap[jAt] = j;
                  jAt++;
               }
            }
            at1IdxMap[iAt] = i;
            iAt++;
         }
      }

//      System.err.printf("\nMol1=%d %d\t mol2=%d %d\n", nAtoms1, atPath1.length,
//               nAtoms2, atPath2.length);
//      printMatrix(distMat);

      int[][] assignment = HungarianAlgorithm.hgAlgorithm(distMat, "min");
      double simSum = 0D;
      for (int i = 0; i < nAtoms1; i++)
      {  double atSim = atSims[at1IdxMap[assignment[i][0]]][at2IdxMap[assignment[i][1]]];
//         System.err.printf("%s\t%s\t%.2g\n",at1IdxMap[assignment[i][0]]+1,at2IdxMap[assignment[i][1]]+1,atSim);
         if( atSim >= minAtSim )
            simSum += atSim;
      }
      return computeMoleculeSimilarity(nAtoms1, nAtoms2, simSum);
   }

   @Override
   double computeAtomSim(int nPaths, int nPaths2, int common)
   {  return common / (double)(Math.max(nPaths, nPaths2) * 2 - common);
   }

   @Override
   double computeMoleculeSimilarity(int nAtoms1, int nAtoms2, double simSum)
   {  return simSum / (2*Math.max(nAtoms1, nAtoms2) - simSum);
   }

   @Override
   public void printMatrix(double[][] m)
   {  for (int i = 0; i < m.length; i++)
      {  double[] row = m[i];
          for (int j = 0; j < row.length; j++)
              System.err.printf("%5.4f ", row[j]);
          System.err.println();
      }
  }

   public static void main(String argv[])
   {  String smiPat;
      String smiTar;

      // aa paper 1
      smiPat = "c1ccn[nH]1";
      smiTar = "Cc1ccno1";

      // aa paper 2
      smiPat = "c1ccnn1C(=O)N1CCC1";
      smiTar = "c1cc[nH]c1C(=O)N1CCC1";
      smiTar = "c1cc[nH]c1";

      smiPat = "c1cccn1C(=O)N1CCC1";
      smiTar = "c1cc[nH]c1";

      OEGraphMol mol1 = new OEGraphMol();
      oechem.OEParseSmiles(mol1, smiPat);
      oechem.OEAssignAromaticFlags(mol1);

      System.err.println(mol1.IsValid() + " " + mol1.NumAtoms()+"\t"+oechem.OECreateAbsSmiString(mol1));

      OEGraphMol mol2 = new OEGraphMol();
      oechem.OEParseSmiles(mol2, smiTar);
      oechem.OEAssignAromaticFlags(mol2);
      System.err.println(mol2.IsValid() + " " + mol2.NumAtoms()+"\t"+oechem.OECreateAbsSmiString(mol2));

      AAPathComparatorFact fact = new AAPathComparatorFact(AAPathCompareType.DEFAULT, 4);
      SimComparator<OEMolBase> simer = fact.createComparator(mol1);
      try{  Thread.sleep(50);
      } catch (InterruptedException e)
      {  e.printStackTrace();
      }
      long start = System.currentTimeMillis();
      double sim = simer.similarity(mol2);
      System.out.printf("%f\n", sim);

      System.exit(0);

      for(int i=0; i<400; i++)
      {  // check if similarity is atom order independent
         oechem.OEScrambleMolecule(mol2);
         double sim2 = simer.similarity(mol2);
         if( Math.abs(sim2 - sim)>0.0001 )
         {  System.out.printf("%f\n", sim2);
            sim=sim2;
         }
      }
      System.err.printf("Done in %.3f\n", (System.currentTimeMillis()-start)/1000D);
      simer.close();
   }
}
