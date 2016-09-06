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

import java.util.Arrays;

import openeye.oechem.*;

import com.aestel.chemistry.openEye.SimComparator;

/**
 * Modified algorithm to compute AtomAtomPath similarity.
 *
 * The algorithm to compute the similarity of two atoms is changed as follows:
 *    pathsInCommon / (2 * max(nPath(AtomMolA), nPath(AtomMolB)) - pathsInCommon)
 *
 * This should prevent the situation where a small substituent results in a higher
 * similarity than a big substituent.
 *
 * This implementation uses int as the type to store the hascode for a path.
 *
 * @author albertgo
 *
 */
public class IAAPathComparatorInt implements SimComparator<OEMolBase>, IAAPathComputerInterface
{  private static final boolean DEBUG = false;

   protected final int maxBonds;
   protected final int[] atomTypes;
   protected final OEAtomBase[] atoms;

   private final OEMolBase mol;
   private final int nAtoms;
   private final int[][] atPath;
   private final HeadAtomComputer headAtomComputer;



   IAAPathComparatorInt(OEMolBase mol, HeadAtomComputer haComputer, int maxBonds)
   {  this.mol = new OEGraphMol(mol);
      this.headAtomComputer = haComputer;
      oechem.OESuppressHydrogens(this.mol);
      oechem.OEAssignAromaticFlags(this.mol);
      oechem.OEPerceiveSymmetry(this.mol);
      this.maxBonds = maxBonds;
      this.nAtoms = this.mol.NumAtoms();
      this.atoms = new OEAtomBase[this.mol.GetMaxAtomIdx()+1];
      OEAtomBaseIter atIt = this.mol.GetAtoms();
      while(atIt.hasNext())
      {  OEAtomBase at = atIt.next();
         atoms[at.GetIdx()] = at;
      }
      atIt.delete();

      IAAPathGeneratorInt apGenerator = new IAAPathGeneratorInt(this.mol, maxBonds);
      this.atPath = getAtomPaths(apGenerator);
      this.atomTypes = apGenerator.getAtomTypes();
   }

   private int[][] getAtomPaths(IAAPathGeneratorInt apGenerator)
   {  int[][] atPaths = new int[atoms.length][];
      int[][] symClassAPath = new int[atoms.length][];

      OEAtomBaseIter aIt = mol.GetAtoms();
      // get path list starting at each atom;
      while( aIt.hasNext() )
      {  OEAtomBase at = aIt.next();
         int aIdx = at.GetIdx();
         int symClass = at.GetSymmetryClass();

         if( symClassAPath[ symClass ] == null )
            symClassAPath[symClass] = apGenerator.getAtomPaths(at);

         atPaths[aIdx] = symClassAPath[symClass];
      }
      aIt.delete();
      return atPaths;
   }

   @Override
   public double similarity(OEMolBase mol2)
   {  IAAPathComparatorInt c = new IAAPathComparatorInt(mol2, headAtomComputer, maxBonds);
      try
      {  return similarity( c );
      } finally
      {  c.close();
      }
   }


   @Override
   public double similarity(SimComparator<OEMolBase> other)
   {  IAAPathComparatorInt o = (IAAPathComparatorInt)other;
      return similarity(this, o);
   }

   /**
    * Default implementation to be improved.
    */
   @Override
   public double similarity(SimComparator<OEMolBase> other, double minSim)
   { return similarity(other);
   }

   protected double similarity(IAAPathComparatorInt m1, IAAPathComparatorInt m2)
   {  // make sure nAtoms1 < nAtoms2
      if( m1.nAtoms > m2.nAtoms )
      {  IAAPathComparatorInt d = m1;
         m1 = m2;
         m2 = d;
      }

      double[][] atSims = getAtomSimilarities(m1, m2);
      int sortIdx[][] = new int[m1.atPath.length][];
      int[] a1IsAssignedTo = new int[m1.atPath.length];
      int[] a2IsAssignedTo = new int[m2.atPath.length];
      Arrays.fill(a1IsAssignedTo, -1);
      Arrays.fill(a2IsAssignedTo, -1);

      for( int a1=0; a1<m1.atPath.length; a1++)
      {  if( m1.atPath[a1] == null ) continue;

         // indexes of atoms in atPath2 in order of similarity to a1
         sortIdx[a1] = sortIdx(atSims[a1]);
      }

      double simSum = 0D;
      for(int nTryAssign=0; nTryAssign<m1.nAtoms; nTryAssign++)
      {  double maxSim = -.1D;
         int a1 = -1;
         int a2 = -1;

         for(int a1Try=0; a1Try<m1.atPath.length; a1Try++)
         {  if( m1.atPath[a1Try] == null ) continue; // not an atom

            // find a1 which was not assigned yet and has highest similarity to a2
            if( a1IsAssignedTo[a1Try] != -1 ) continue;

            double sim = -1D;
            for(int a2Rank=0; a2Rank<m2.atPath.length; a2Rank++)
            {  int a2Try = sortIdx[a1Try][a2Rank];
               if( a2IsAssignedTo[a2Try] >= 0 ) continue;

               sim = atSims[a1Try][a2Try];

               if( sim > maxSim )
               {  maxSim = sim;
                  a1 = a1Try;
                  a2 = a2Try;
                  break;
               }
            }
         }

         assert a1IsAssignedTo[a1] == -1;
         assert a2IsAssignedTo[a2] == -1;
         a1IsAssignedTo[a1] = a2;
         a2IsAssignedTo[a2] = a1;

         simSum += maxSim;
      }

      return computeMoleculeSimilarity(m1.nAtoms, m2.nAtoms, simSum);
   }


   protected double[][] getAtomSimilarities(IAAPathComparatorInt m1, IAAPathComparatorInt m2)
   {  int[][] atPath1 = m1.atPath;
      int[][] atPath2 = m2.atPath;

      double[][] simMatrix = new double[atPath1.length][atPath2.length];
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
               sim = atomSimilarity(a1, m1, a2, m2);
            simMatrix[a1][a2] = sim;
         }
      }

      if( DEBUG ) printMatrix(simMatrix);
      return simMatrix;
   }

   private double atomSimilarity(int atIdx1, IAAPathComparatorInt m1, int atIdx2, IAAPathComparatorInt m2)
   {  double headAtomSym = headAtomComputer.getHeadAtomSim(atIdx1, m1, atIdx2, m2);
      if( headAtomSym == 0D ) return 0D;

      int[] atPath1 = m1.atPath[atIdx1];
      int[] atPath2 = m2.atPath[atIdx2];

      // take care of disconnected atom that has no paths
      if( atPath1.length == 0 )
      {  if( atPath2.length == 0 ) return headAtomSym;
         return headAtomSym * headAtomSym/(headAtomSym+atPath2.length);
      } else if( atPath2.length == 0 )
      {  return headAtomSym * headAtomSym/(headAtomSym+atPath1.length);
      }

      int m1Features = atPath1.length;
      int m2Features = atPath2.length;
      int common = 0;
      int pIdx = 0;
      int pIdx2 = 0;
      int p = atPath1[pIdx];
      int p2 = atPath2[pIdx2];


      // compute features in common
      // assume feature lists( atPAth1 and atPath2 are sorted
      do
      {  if( p == p2 )
         {  common++;
            if( ++pIdx  >= m1Features || ++pIdx2 >= m2Features ) break;
            p = atPath1[pIdx];
            p2 = atPath2[pIdx2];

         }else if( p < p2 )
         {  if( ++pIdx >= m1Features ) break;
            p = atPath1[pIdx];
         }else // p > p2
         {  if( ++pIdx2 >= m2Features ) break;
            p2 = atPath2[pIdx2];
         }
      }while(true);

      return headAtomSym * computeAtomSim(m1Features, m2Features, common, headAtomSym);
   }


   double computeAtomSim(int nPaths, int nPaths2, int common, double headAtomSym)
   {  // headAtomSim is added to account for similarity of the mapped atoms themselfs
      return (common + headAtomSym) / (Math.max(nPaths, nPaths2) * 2 - common + headAtomSym);
   }


   double computeMoleculeSimilarity(int nAtoms1, int nAtoms2, double simSum)
   {  return simSum / (2*Math.max(nAtoms1, nAtoms2) - simSum);
   }

   @Override
   public int getAtomType(int aIdx)
   {  return atomTypes[aIdx];
   }

   @Override
   public OEAtomBase getAtom(int aIdx)
   {  return atoms[aIdx];
   }

   @Override
   public int getAtomNum(int atomType)
   {  return IAAPathGenerator.atomTypeToAtomNum(atomType);
   }

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
   private static boolean less(final double x, final double y)
   {  return x > y;
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

   @Override
   public void close()
   {  mol.delete();
   }
}
