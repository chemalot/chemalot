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

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import openeye.oechem.*;

import com.aestel.chemistry.openEye.SimComparator;
import com.aestel.chemistry.openEye.tools.DFSIterator;
import com.aestel.chemistry.openEye.tools.OEAtomBondPath;
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
 * The total similarity is computed as
 *    sum(Matched(atomPairSimilarity))/ (nAtom1 + nAtom2 - sum(M(aps))
 *
 * @author albertgo 2010
 *
 */
public class AAPathComparator implements SimComparator<OEMolBase>
{  private static final boolean DEBUG = false;

   private final OEMolBase mol;
   final int nAtoms;
   final String[][] atPath;

   protected final double minAtSim;
   private final int maxBonds;

   /**
    *
    * @param mol used to compute similarity
    * @param maxBonds number of bonds for linear paths
    * @param minAtSim minimum similarity of other atom to be considered
    */
   AAPathComparator(OEMolBase mol, int maxBonds, double minAtSim)
   {  this.mol = new OEGraphMol(mol);
      oechem.OESuppressHydrogens(this.mol);
      oechem.OEPerceiveSymmetry(this.mol);
      this.nAtoms = this.mol.NumAtoms();

      this.minAtSim = minAtSim;
      this.maxBonds = maxBonds;
      this.atPath = getAtomPaths(this.mol);
   }

   /** Get the pathsets originating from each atom in lMol.
    *
    * Specialized implementation is 20% faster but is it worth it?
    * This uses a hand crafted depth first search instead of the {@see DFSIterator}
    *
    */
   private String[][] getAtomPaths(OEMolBase lMol)
   {  AtomPathGenerator apGenerator = new AtomPathGenerator(lMol, maxBonds);
      String[][] atPaths = new String[lMol.GetMaxAtomIdx()+1][];
      String[][] symClassAPath = new String[lMol.GetMaxAtomIdx()+1][];
      OEAtomBaseIter aIt = lMol.GetAtoms();

      // get path list starting at each atom;
      while( aIt.hasNext() )
      {  OEAtomBase at = aIt.next();
         int aIdx = at.GetIdx();

         if( symClassAPath[ aIdx ] == null )
         {  String[] singleAtPaths = apGenerator.getAtomPaths(at);
            Arrays.sort(singleAtPaths);
            symClassAPath[aIdx] = singleAtPaths;
         }
         atPaths[aIdx] = symClassAPath[aIdx];
      }
      aIt.delete();
      return atPaths;
   }

   /* Gives same results but is 20% slower.
   * This uses a well tested DFS implementation.
   */
  @SuppressWarnings("unused")
   private String[][] getAtomPathsO(OEMolBase lMol)
   {  String[][] atPaths = new String[lMol.GetMaxAtomIdx()+1][];
      OEAtomBaseIter aIt = lMol.GetAtoms();
      ArrayList<String> pathList = new ArrayList<String>();

      // get path list starting at each atom;
      while( aIt.hasNext() )
      {  OEAtomBase at = aIt.next();
         DFSIterator dfIt = new DFSIterator(at, maxBonds);
         while(dfIt.hasNext())
         {  OEAtomBondPath path = dfIt.next();
            pathList.add(path.toString());
         }
         dfIt.close();

         String[] singleAtPaths = pathList.toArray(new String[pathList.size()]);
         Arrays.sort(singleAtPaths);
         atPaths[at.GetIdx()] = singleAtPaths;
         pathList.clear();
      }
      aIt.delete();
      return atPaths;
   }

   @Override
   public double similarity(OEMolBase mol2)
   {  OEGraphMol tmp = new OEGraphMol(mol2);
      oechem.OESuppressHydrogens(tmp);
      String[][] atPath2 = getAtomPaths(tmp);
      int nAtoms2 = tmp.NumAtoms();
      tmp.delete();
      return similarity( nAtoms, atPath, nAtoms2, atPath2 );
   }


   /**
    * Conpute similarity of two molecules given pathsets.
    * @param nAtoms1 number of atoms in 1
    * @param atPath1 pathsets for 1
    * @param nAtoms2 number of atoms in 2
    * @param atPath2 pathsets for 1
    */
   protected double similarity(int nAtoms1, String[][] atPath1,
                             int nAtoms2, String[][] atPath2)
   {  // make user atPath length is smaller than atPath2 length,
      if( atPath1.length > atPath2.length )
      {  // swap
         String[][] d = atPath1;
         atPath1 = atPath2;
         atPath2= d;
      }

      double[][] atSims = getAtomSimilarities(atPath1, atPath2);
      int sortIdx[][] = new int[atPath1.length][];
      int[] a1IsAssignedTo = new int[atPath1.length];
      int[] a2IsAssignedTo = new int[atPath2.length];
      Arrays.fill(a1IsAssignedTo, -1);
      Arrays.fill(a2IsAssignedTo, -1);

      for( int a1=0; a1<atPath1.length; a1++)
      {  if( atPath1[a1] == null ) continue;

         // indexes of atoms in atPath2 in order of similarity to a1
         sortIdx[a1] = sortIdx(atSims[a1]);
      }

      double simSum = 0D;
      for(int nTryAssign=0; nTryAssign<nAtoms1; nTryAssign++)
      {  double maxSim = -.1D;
         int a1 = -1;
         int a2 = -1;

         for(int a1Try=0; a1Try<atPath1.length; a1Try++)
         {  if( atPath1[a1Try] == null ) continue; // this is an empty oepeneye atom

            if( a1IsAssignedTo[a1Try] != -1 ) continue; // already assigned

            // find a1 which was not assigned yet and has highest similarity to a2
            double sim = -1D;
            for(int a2Rank=0; a2Rank<atPath2.length; a2Rank++)
            {  int a2Try = sortIdx[a1Try][a2Rank];
               if( a2IsAssignedTo[a2Try] >= 0 ) continue;

               sim = atSims[a1Try][a2Try];

               if( sim < minAtSim ) break; // no atoms are similar enough

               if( sim > maxSim )
               {  maxSim = sim;
                  a1 = a1Try;
                  a2 = a2Try;
                  break;
               }
            }
         }

         if( maxSim >= minAtSim) // atom assignment accomplished: similar enough
         {  assert a1IsAssignedTo[a1] == -1;
            assert a2IsAssignedTo[a2] == -1;
            a1IsAssignedTo[a1] = a2;
            a2IsAssignedTo[a2] = a1;

            simSum += maxSim;
         }
      }

      return computeMoleculeSimilarity(nAtoms1, nAtoms2, simSum);
   }

   double computeMoleculeSimilarity(int nAtoms1, int nAtoms2, double simSum)
   {  return simSum / (nAtoms1 + nAtoms2 - simSum);
   }


   protected double[][] getAtomSimilarities(String[][] atPath1, String[][] atPath2)
   {  double[][] simMatrix = new double[atPath1.length][atPath2.length];

      for(int a1=0; a1<atPath1.length; a1++)
      {  if( atPath1[a1] == null )   // oechem might have atom indexes which do not exist
         {  Arrays.fill(simMatrix[a1], -.1D);
            continue;
         }

         for(int a2=0; a2<atPath2.length; a2++)
         {  double sim;
            if( atPath2[a2] == null )
               sim = -.1D;
            else
               sim = atomSimilarity(atPath1[a1], atPath2[a2]);
            simMatrix[a1][a2] = sim;
         }
      }

      if( DEBUG ) printMatrix(simMatrix);
      return simMatrix;
   }

   private double atomSimilarity(String[] paths, String[] paths2)
   {  int common = 0;
      int pIdx = 0;
      int pIdx2 = 0;
      String p = paths[pIdx];
      String p2 = paths2[pIdx2];

      // assume path's are sorted alphabetically
      // go through list and:
      //  increment mismatch counter on missmatch
      //  increment pIdx or pIdx depending on which string is smaller
      do
      {  int compare = p.compareTo(p2);
         if( compare == 0 )
         {  common++;
            if( ++pIdx  >= paths.length  ) break;
            if( ++pIdx2 >= paths2.length ) break;
            p = paths[pIdx];
            p2 = paths2[pIdx2];

         }else if( compare < 0 )
         {  if( ++pIdx >= paths.length ) break;
            p = paths[pIdx];
         }else
         {  if( ++pIdx2 >= paths2.length ) break;
            p2 = paths2[pIdx2];
         }
      }while(true);

      return computeAtomSim(paths.length, paths2.length, common);
   }


   double computeAtomSim(int nPaths, int nPaths2, int common)
   {  return common / (double)(nPaths + nPaths2 - common);
   }

   @Override
   public double similarity(SimComparator<OEMolBase> other)
   {  AAPathComparator o = (AAPathComparator)other;
      return similarity(this.nAtoms, this.atPath, o.nAtoms, o.atPath);
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
   {  mol.delete();
   }

   /**
    * return an array of indexes into ds such that they point to an decreasing
    * sequence of values.
    */
   private static int[] sortIdx(double[] ds)
   {  int[] indexes = new int[ds.length];
      for(int i=0; i<indexes.length; i++)
         indexes[i] = i;

      double[] tmp = ds.clone(); // qs changes ds
      quicksort(tmp, indexes, 0, indexes.length - 1);

      return indexes;
   }

   // quicksort a[left] to a[right]
   private static void quicksort(double[] a, int[] index, int left, int right)
   {  if (right <= left)
         return;
      int i = partition(a, index, left, right);
      quicksort(a, index, left, i - 1);
      quicksort(a, index, i + 1, right);
   }

   // partition a[left] to a[right], assumes left < right
   private static int partition(double[] a, int[] index, int left, int right)
   {  int i = left - 1;
      int j = right;
      while (true)
      {
         while (less(a[++i], a[right]))
            // find item on left to swap
            ; // a[right] acts as sentinel
         while (less(a[right], a[--j]))
            // find item on right to swap
            if (j == left)
               break; // don't go out-of-bounds
         if (i >= j)
            break; // check if pointers cross
         exch(a, index, i, j); // swap two elements into place
      }
      exch(a, index, i, right); // swap with partition element
      return i;
   }

   // is x < y ?
   private static boolean less(double x, double y)
   {  return (x > y);
   }

   // exchange a[i] and a[j]
   private static void exch(double[] a, int[] index, int i, int j)
   {  double swap = a[i];
      a[i] = a[j];
      a[j] = swap;
      int b = index[i];
      index[i] = index[j];
      index[j] = b;
   }

   public void printMatrix(double[][] m)
   {  for (int i = 0; i < m.length; i++)
      {  double[] row = m[i];
          for (int j = 0; j < row.length; j++)
              System.err.printf("%5.4f ", row[j]);
          System.err.println();
      }
  }

   public static void main(String argv[]) throws IOException
   {  String smiPat;
      String smiTar;
      String smiTar2;

      // aa paper 1
      smiPat = "c1ccn[nH]1";
      smiTar = "Cc1ccno1";

      // aa paper 2
      smiPat = "c1ccnn1C(=O)N1CCC1";
      smiTar = "c1cc[nH]c1C(=O)N1CCC1"; // 0.321
      smiTar2 = "c1ccn[nH]1";       // 0.193
      smiTar = "c1cc[nH]c1";       // 0.057

      smiPat = "c1cccn1C(=O)N1CCC1";
      smiTar = "c1cc[nH]c1";       // 0.194

      OEGraphMol mol1 = new OEGraphMol();
      oechem.OEParseSmiles(mol1, smiPat);
      oechem.OEAssignAromaticFlags(mol1);

      System.err.println(mol1.IsValid() + " " + mol1.NumAtoms()+"\t"+oechem.OECreateAbsSmiString(mol1));

      OEGraphMol mol2 = new OEGraphMol();
      oechem.OEParseSmiles(mol2, smiTar);
      oechem.OEAssignAromaticFlags(mol2);
      System.err.println(mol2.IsValid() + " " + mol2.NumAtoms()+"\t"+oechem.OECreateAbsSmiString(mol2));

      OEGraphMol mol3 = new OEGraphMol();
      oechem.OEParseSmiles(mol3, smiTar2);
      oechem.OEAssignAromaticFlags(mol3);

      AAPathComparatorFact fact = new AAPathComparatorFact(AAPathCompareType.FUZZY, AAPathComparatorFact.DEFAULTVersion);
//      AAPathComparatorFact fact = new AAPathComparatorFact(AAPathCompareType.DEFAULT, AAPathComparatorFact.DEFAULTVersion);
      SimComparator<OEMolBase> simer = fact.createComparator(mol1);
      try{  Thread.sleep(50);
      } catch (InterruptedException e)
      {  e.printStackTrace();
      }
      double sim = simer.similarity(mol2);
      System.out.printf("%f\n", sim);
      sim = simer.similarity(mol3);
      System.out.printf("%f\n", sim);

      //System.exit(0);
      //System.err.print("Waiting, hit enter: ");System.in.read();
      long start = System.currentTimeMillis();
      SimComparator<OEMolBase> compa2 = fact.createComparator(mol2);
      for(int i=0; i<20000; i++)
      {  // check if similarity is atom order independent
         //oechem.OEScrambleMolecule(mol2);
         // slower atom paths are computed:
         //double sim2 = simer.similarity(mol2);

         // faster: atompaths are pre computed
         double sim2 = simer.similarity(compa2);

         if( Math.abs(sim2 - sim)>0.0001 )
         {  System.out.printf("%f\n", sim2);
            sim=sim2;
         }
      }
      System.err.printf("Done in %.3f\n", (System.currentTimeMillis()-start)/1000D);
      simer.close();
   }
}
