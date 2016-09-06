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

import java.util.ArrayList;
import java.util.Arrays;

import openeye.oechem.*;

import com.aestel.chemistry.openEye.SimComparator;
import com.genentech.oechem.tools.OETools;

/**
 * Modified algorithm to compute AtomAtomPath similarity.
 *
 * The algorithm to compute the similarity of two atoms is changed as follows:
 *    pathsInCommon / (2 * max(nPath(AtomMolA), nPath(AtomMolB)) - pathsInCommon)
 *
 * This should prevent the situation where a small substituent results in a higher
 * similarity than a big substituent.
 *
 * This implementation uses char as the type to store the hascode for a path.
 *
 * @author albertgo
 *
 */
public class IAAPathComparatorChar implements SimComparator<OEMolBase>, IAAPathComputerInterface
{  private static final boolean DEBUG = false;

   protected final int maxBonds;
   protected final int[] atomTypes;
   protected final OEAtomBase[] atoms;

   private final OEMolBase mol;
   private final int nAtoms;
   private final char[][] atPath;
   private final HeadAtomComputer headAtomComputer;



   IAAPathComparatorChar(OEMolBase mol, HeadAtomComputer haComputer, int maxBonds)
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

      IAAPathGeneratorChar apGenerator = new IAAPathGeneratorChar(this.mol, maxBonds);
      this.atPath = getAtomPaths(apGenerator);
      this.atomTypes = apGenerator.getAtomTypes();
   }

   private char[][] getAtomPaths(IAAPathGeneratorChar apGenerator)
   {  char[][] atPaths = new char[atoms.length][];
      char[][] symClassAPath = new char[atoms.length][];

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
   {  IAAPathComparatorChar c = new IAAPathComparatorChar(mol2, headAtomComputer, maxBonds);
      try
      {  return similarity( c );
      } finally
      {  c.close();
      }
   }


   @Override
   public double similarity(SimComparator<OEMolBase> other)
   {  IAAPathComparatorChar o = (IAAPathComparatorChar)other;
      return similarity(this, o);
   }


   /**
    * Default implementation to be improved.
    */
   @Override
   public double similarity(SimComparator<OEMolBase> other, double minSim)
   {  IAAPathComparatorChar o = (IAAPathComparatorChar)other;
      if((Math.min(this.nAtoms,o.nAtoms)*10000)/Math.max(this.nAtoms,o.nAtoms)<(int)(minSim*10000D))
      {  return 0D;
      }
      return similarity(other);
   }

   /* idea to improve mapping for tied atoms:
    * discussed on 1/2/2015 at Mandalay (ML,AG)
    * 1) find highest similarity atom pair (a,b) -> map, mark as visited
    * 2) put neighbors of a into nighbor list A, and neighbors of b into neighbor list B
    * 3) look for highest similarity pair (a1,b1) between atoms in A and B
    * 4) if sim(a1,b1) == 0: mark members of A as visited, goto 1
    * 5) map a1 -> b1
    * 6) remove a1 from A and b1 from B
    * 7) if not all atoms visited goto 2
    * 8) check visited but not mapped atoms for mappings
    *
    * Variant:
    *   - in case of ties try each tie and keep highest sim mapping.
    *   - in step 4 use "< 0.05" instead of "== 0"
    *
    */




   protected double similarity(IAAPathComparatorChar m1, IAAPathComparatorChar m2)
   {  // make sure nAtoms1 < nAtoms2
      if( m1.nAtoms > m2.nAtoms )
      {  IAAPathComparatorChar d = m1;
         m1 = m2;
         m2 = d;
      }

      // atPairs are sorted by similarity
      AtomPair[] atPairs = getAtomSimilarities(m1, m2);
      if( DEBUG )
      {  System.err.printf("%s\t%s\n",
               OETools.molToCanSmi(m1.mol,true), OETools.molToCanSmi(m2.mol,true));
         for( AtomPair ap : atPairs ) System.err.println(ap);
      }
      int[] a1IsAssignedTo = new int[m1.atPath.length];
      int[] a2IsAssignedTo = new int[m2.atPath.length];
      Arrays.fill(a1IsAssignedTo, -1);
      Arrays.fill(a2IsAssignedTo, -1);

      double simSum = 0D;
      for(int mostSimilarIdx = 0; mostSimilarIdx < atPairs.length; mostSimilarIdx++)
      {  int a1 = atPairs[mostSimilarIdx].at1Idx;
         if( a1IsAssignedTo[a1] != -1 ) continue;

         int a2 = atPairs[mostSimilarIdx].at2Idx;
         if( a2IsAssignedTo[a2] != -1 ) continue;

         double sim = atPairs[mostSimilarIdx].sim;
         simSum += sim;

         a1IsAssignedTo[a1] = a2;
         a2IsAssignedTo[a2] = a1;
      }

      return computeMoleculeSimilarity(m1.nAtoms, m2.nAtoms, simSum);
   }


   /** @return array of atom pairs with similarity > 0 */
   protected AtomPair[] getAtomSimilarities(IAAPathComparatorChar m1, IAAPathComparatorChar m2)
   {  char[][] atPath1 = m1.atPath;
      char[][] atPath2 = m2.atPath;

      ArrayList<AtomPair> atPairs = new ArrayList<AtomPair>(atPath1.length * atPath2.length);

      for(int a1=0; a1<atPath1.length; a1++)
      {  if( atPath1[a1] == null )   // oechem might have atom indexes which do not exist
            continue;

         for(int a2=0; a2<atPath2.length; a2++)
         {  double sim;
            if( atPath2[a2] != null )
            {  sim = atomSimilarity(a1, m1, a2, m2);
               if( sim > 0D )
                  atPairs.add(new AtomPair(a1, a2, sim));
            }
         }
      }


      AtomPair[] ret = atPairs.toArray(new AtomPair[atPairs.size()]);
      Arrays.sort(ret);
      return ret;
   }

   private double atomSimilarity(int atIdx1, IAAPathComparatorChar m1, int atIdx2, IAAPathComparatorChar m2)
   {  double headAtomSim = headAtomComputer.getHeadAtomSim(atIdx1, m1, atIdx2, m2);
      if( headAtomSim == 0D ) return 0D;

      char[] atPath1 = m1.atPath[atIdx1];
      char[] atPath2 = m2.atPath[atIdx2];

      // take care of disconnected atom that has no paths
      if( atPath1.length == 0 )
      {  if( atPath2.length == 0 ) return headAtomSim;
         return headAtomSim * headAtomSim/(headAtomSim+atPath2.length);
      } else if( atPath2.length == 0 )
      {  return headAtomSim * headAtomSim/(headAtomSim+atPath1.length);
      }

      int m1Features = atPath1.length;
      int m2Features = atPath2.length;
      int common = 0;
      int pIdx = 0;
      int pIdx2 = 0;
      char p = atPath1[pIdx];
      char p2 = atPath2[pIdx2];


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

      return headAtomSim * computeAtomSim(m1Features, m2Features, common, headAtomSim);
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


   @Override
   public void close()
   {  mol.delete();
   }
}

class AtomPair implements Comparable<AtomPair>
{  AtomPair(int at1Idx, int at2Idx, double sim)
   {  this.at1Idx = at1Idx;
      this.at2Idx = at2Idx;
      this.sim = sim;
   }

   int at1Idx;
   int at2Idx;
   double sim;


   @Override
   public int compareTo(AtomPair o)
   {  if( this.sim > o.sim ) return -1;
      if( this.sim < o.sim ) return 1;
      if( this.at2Idx < o.at2Idx ) return -1;
      if( this.at2Idx > o.at2Idx ) return 1;
      if( this.at1Idx < o.at1Idx ) return -1;
      if( this.at1Idx > o.at1Idx ) return 1;
      return 0;
   }

   @Override
   public String toString()
   {  return String.format("%2d %2d %.4f", at1Idx, at2Idx, sim);
   }
}
