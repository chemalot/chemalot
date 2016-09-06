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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.TimeoutException;

import openeye.oechem.*;

import com.aestel.chemistry.molecule.Atom;
import com.aestel.chemistry.openEye.tools.DFSIterator;
import com.aestel.chemistry.openEye.tools.OEAtomBondPath;



/**
 * Compute topological indexes.
 *
 * References:
 * (1) Balaban, A. T. Highly discriminating distance-based topological index. Chemical Physics Letters 1982, 89, 399-404.
 * (2) Balaban, A. T. Topological indices based on topological distances in molecular graphs. Pure and Applied Chemistry 1983, 55, 199-206.
 * (3) Bonchev, D.; Mountain, C. F.; Seitz, W. A.; Balaban, A. T. Modeling the anticancer action of some retinoid compounds by making use of the OASIS method. J. Med. Chem. 1993, 36, 1562-1569.
 * (4) Ivanciuc, O.; Ivanciuc, T.; Balaban, A. T. Design of Topological Indices. Part 10.1 Parameters Based on Electronegativity and Covalent Radius for the Computation of Molecular Graph Descriptors for Heteroatom-Containing Molecules. J. Chem. Inf. Comput. Sci. 1998, 38, 395-401.
 * (5) Todeschini, R.; Consonni, V. Handbook of Molecular Descriptors 2008, 23
 * (6) Mercader, A.; Castro, E. A.; Toropov, A. A. Maximum Topological Distances Based Indices as Molecular Descriptors for QSPR. 4. Modeling the Enthalpy of Formation of Hydrocarbons from Elements. International Journal of Molecular Sciences 2001, 2, 121-132.
 * @author albertgo
 *
 */

public class TopologicalIndexer
{  private OEGraphMol mol = new OEGraphMol();
   private int nAtom;
   private int nBond;
   private OEAtomBase[] ats;
   private int[][] dMat;
   private int[][] multiDMat12;
   private int[][] multiLongDMat12;
   private long[] distanceDegrees;
   private int[] vertexDegrees;
   private int[] distributionDistDegrees;
   private int[][] adjacentPairs;
   private double[] multiDegrees;
   private double[] multiLongDegrees;


   public TopologicalIndexer(  )
   {
   }

   /**
    * call this method before any of the get>>> functions
    * @param inmol
    */
   public void computeIndexes(OEMolBase inmol)
   {  mol.Clear();
      oechem.OEAddMols(mol, inmol);
      oechem.OESuppressHydrogens(mol);
      oechem.OEAssignAromaticFlags(mol);

      nAtom = mol.NumAtoms();
      nBond = mol.NumBonds();
      ats  = null;
      dMat = null;
      adjacentPairs = null;
      vertexDegrees = null;
      distanceDegrees = null;
      distributionDistDegrees = null;
      multiDMat12  = null;
      multiLongDMat12 = null;
      multiDegrees  = null;
      multiLongDegrees  = null;
   }


   public double getBalabanJIndex() throws TimeoutException
   {  getAdjacentPairs();
      getMultiDegrees();
      double ret = 0D;

      for(int i=0; i<adjacentPairs.length; i++)
      {  int[] pair = adjacentPairs[i];
         double s1 = multiDegrees[pair[0]];
         double s2 = multiDegrees[pair[1]];
         ret += 1/Math.sqrt(s1*s2);
      }

      return (double)nBond/(nBond-nAtom + 1 +1) * ret;
   }


   /** Eq. 6 and 8 in Ref 3
    * @throws TimeoutException if a single molecule needs more thatn 30sec */
   public double getBalabanJYIndex() throws TimeoutException
   {  double ret = 0D;
      getAdjacentPairs();
      getMultiDegrees();

      for(int i=0; i<adjacentPairs.length; i++)
      {  int[] pair = adjacentPairs[i];
         double s1 = getBalabanDY(pair[0], multiDegrees[pair[0]]);
         double s2 = getBalabanDY(pair[1], multiDegrees[pair[1]]);
         ret += 1/Math.sqrt(s1*s2);
      }

      return (double)nBond/(nBond-nAtom + 1 +1) * ret;
   }

   /** Eq. 6 and 8 in Ref 3
    * @throws TimeoutException if a single molecule needs more thatn 30sec */
   public double getBalabanJXIndex() throws TimeoutException
   {  double ret = 0D;
      getAdjacentPairs();
      getMultiDegrees();

      for(int i=0; i<adjacentPairs.length; i++)
      {  int[] pair = adjacentPairs[i];
         double s1 = getBalabanDX(pair[0], multiDegrees[pair[0]]);
         double s2 = getBalabanDX(pair[1], multiDegrees[pair[1]]);
//System.err.println(ret + " " + s1 + " " + s2+ " " + Math.sqrt(s1*s2));
         ret += 1/Math.sqrt(s1*s2);
      }

      return (double)nBond/(nBond-nAtom + 1 +1) * ret;
   }


   public double getBalabanJStarIndex() throws TimeoutException
   {  getAdjacentPairs();
      getMultiLongDegrees();
      double ret = 0D;

      for(int i=0; i<adjacentPairs.length; i++)
      {  int[] pair = adjacentPairs[i];
         double s1 = multiLongDegrees[pair[0]];
         double s2 = multiLongDegrees[pair[1]];
         ret += 1/Math.sqrt(s1*s2);
      }

      return (double)nBond/(nBond-nAtom + 1 +1) * ret;
   }


   /** Eq. 6 and 8 in Ref 3
    * @throws TimeoutException if a single molecule needs more thatn 30sec */
   public double getBalabanJYStarIndex() throws TimeoutException
   {  getAdjacentPairs();
      getMultiLongDegrees();
      double ret = 0D;

      for(int i=0; i<adjacentPairs.length; i++)
      {  int[] pair = adjacentPairs[i];
         double s1 = getBalabanDY(pair[0], multiLongDegrees[pair[0]]);
         double s2 = getBalabanDY(pair[1], multiLongDegrees[pair[1]]);
         ret += 1/Math.sqrt(s1*s2);
      }

      return (double)nBond/(nBond-nAtom + 1 +1) * ret;
   }

   /** Eq. 6 and 8 in Ref 3
    * @throws TimeoutException if a single molecule needs more thatn 30sec */
   public double getBalabanJXStarIndex() throws TimeoutException
   {  getAdjacentPairs();
      getMultiLongDegrees();
      double ret = 0D;

      for(int i=0; i<adjacentPairs.length; i++)
      {  int[] pair = adjacentPairs[i];
         double s1 = getBalabanDX(pair[0], multiLongDegrees[pair[0]]);
         double s2 = getBalabanDX(pair[1], multiLongDegrees[pair[1]]);
         ret += 1/Math.sqrt(s1*s2);
      }

      return (double)nBond/(nBond-nAtom + 1 +1) * ret;
   }


   public long getWienerIndex()
   {  long ret = 0;
      getDistanceDegree();

      for(int i=0; i<distanceDegrees.length; i++)
         ret += distanceDegrees[i];

      return ret/2;
   }

   /** M1 as in ref 3 */
   public int getZagrebIndex()
   {  int ret = 0;
      getVertexDegrees();

      for(int i=0; i<vertexDegrees.length; i++)
         ret += vertexDegrees[i] * vertexDegrees[i];

      return ret;
   }


   public void close()
   {  mol.delete();
   }


   @SuppressWarnings({ "rawtypes", "unchecked" })
   private int[][] getAdjacentPairs()
   {  if( adjacentPairs != null ) return adjacentPairs;
      getDistanceMatrix();

      List pairs = new ArrayList();
      for(int i=0; i<dMat.length; i++)
         for(int j=0; j<i; j++)
            if(dMat[i][j] == 1) pairs.add(new int[] {i,j});

      adjacentPairs = (int[][]) pairs.toArray(new int[pairs.size()][2]);
      return adjacentPairs;
   }


   @SuppressWarnings("unused")
   private int[] getDistribution()
   {  if( distributionDistDegrees != null ) return distributionDistDegrees;
      getDistanceMatrix();

      int[] dist = new int[dMat.length];
      for(int d=1; d < distributionDistDegrees.length; d++)
         for(int i=0; i<dMat.length; i++)
            for(int j=0; j<i; j++)
               if( dMat[i][j] == d)
                  dist[d]++;

      int maxDeg = 0;
      for(int i=1; i< dist.length; i++)
         if(dist[i] > 0 ) maxDeg = i;

      distributionDistDegrees = new int[maxDeg+1];
      System.arraycopy(dist, 0, distributionDistDegrees, 0, maxDeg+1);
      return distributionDistDegrees;
   }

   /** array of number of bonds per atom by atom index */
   private int[] getVertexDegrees()
   {  if( vertexDegrees != null ) return vertexDegrees;

      vertexDegrees = new int[mol.GetMaxAtomIdx()+1];
      OEAtomBaseIter atIt = mol.GetAtoms();
      while(atIt.hasNext())
      {  OEAtomBase at = atIt.next();
         vertexDegrees[at.GetIdx()] = at.GetExplicitDegree();
      }
      atIt.delete();

      return vertexDegrees;
   }


   private long [] getDistanceDegree()
   {  if( distanceDegrees != null ) return distanceDegrees;
      getDistanceMatrix();

      distanceDegrees = getDegrees(dMat);

      return distanceDegrees;
   }


   /**
    * MultiDegree of an atom is the sum e of its row (=column) in the multi distance matrix.
    * @throws TimeoutException if a single molecule needs more thatn 30sec
    */
   private double[] getMultiDegrees() throws TimeoutException
   {  if( multiDegrees != null ) return multiDegrees;
      getMultiDistanceMatrixTimes12();

      long[] deg = getDegrees(multiDMat12);
      multiDegrees = new double[deg.length];
      for( int i=0; i<deg.length; i++)
         multiDegrees[i] = deg[i]/(double)12;

      return multiDegrees;
   }

   /**
    * MultiDegree of an atom is the sum e of its row (=column) in the multi distance matrix.
    * @throws TimeoutException if a single molecule needs more thatn 30sec
    */
   private double[] getMultiLongDegrees() throws TimeoutException
   {  if( multiLongDegrees != null ) return multiLongDegrees;
      getMultiDistanceMatrixTimes12();

      long[] deg = getDegrees(multiLongDMat12);
      multiLongDegrees = new double[deg.length];
      for( int i=0; i<deg.length; i++)
         multiLongDegrees[i] = deg[i]/(double)12;

      return multiLongDegrees;
   }


   /**
    * degree of an atom is the sum of its row (=column) in the distance matrix.
    */
   private long[] getDegrees(int[][] dMat2)
   {  getAtoms();
      long[] deg = new long[dMat2.length];
      for(int i=0; i<dMat2.length; i++)
      {  if(ats[i] == null) continue;
         for(int j=0; j<dMat2.length; j++)
         {  if(ats[j] == null) continue;
            deg[i] += dMat2[i][j];
         }
      }
      return deg;
   }


   /**
    * fill ats array with atoms indexed by at.getIdx()
    *
    */
   OEAtomBase[] getAtoms()
   {  if( ats != null ) return ats;

      ats = new OEAtomBase[mol.GetMaxAtomIdx() + 1];
      OEAtomBaseIter atIt1 = mol.GetAtoms();
      while (atIt1.hasNext())
      {
         OEAtomBase at = atIt1.next();
         ats[at.GetIdx()] = at;
      }
      atIt1.delete();

      return ats;
   }

   /**
    * Calculate distance matrix of atoms in molecule indexed by {@see OEAtomBase#getIndex}.
    */
   private int[][] getDistanceMatrix()
   {  if( dMat != null ) return dMat;
      getAtoms();

      dMat = new int[ats.length][ats.length];
      for( int i=0; i<ats.length; i++)
      {  if( ats[i] == null) continue;
         for(int j=0; j<i; j++)
         {  if( ats[j] == null ) continue;
            dMat[i][j] = oechem.OEGetPathLength(ats[i], ats[j]);
            dMat[j][i] = dMat[i][j];
         }
      }

      return dMat;
   }


   /**
    * Fill multiDMat12 with the MultiBondDistance Matrix times 12 (cf. reference 1).
    *
    *  And the multiShortDMat12 which is called delta in ref. 6
    * @throws TimeoutException if a single molecule needs more thatn 30sec. Eg. for buckyballs.
    */
   private int[][] getMultiDistanceMatrixTimes12() throws TimeoutException
   {  int maxAtIdx = mol.GetMaxAtomIdx();
      int nBonds = mol.NumBonds();
      int[] multiBondOrdersTimes12 = new int[mol.GetMaxBondIdx()+1];

      OEBondBaseIter bdIt = mol.GetBonds();
      while(bdIt.hasNext())
      {  OEBondBase bd = bdIt.next();
         multiBondOrdersTimes12[bd.GetIdx()] = getMultiBondOrderTimes12(bd);
      }
      bdIt.delete();

      multiDMat12 = new int[maxAtIdx+1][maxAtIdx+1];
      multiLongDMat12 = new int[maxAtIdx+1][maxAtIdx+1];

      OEAtomBaseIter atIt = mol.GetAtoms();
      while(atIt.hasNext())
      {  OEAtomBase at = atIt.next();
         int i = at.GetIdx();
         Arrays.fill(multiDMat12[i], Integer.MAX_VALUE);
         Arrays.fill(multiLongDMat12[i], Integer.MIN_VALUE);

         DFSIterator dfs = new DFSIterator(at, nBonds);
         long start = System.currentTimeMillis();
         while( dfs.hasNext())
         {  if(System.currentTimeMillis() - start > 30*1000)
            {  dfs.close();
               atIt.delete();
               throw new TimeoutException("Exceeded time in topological Indexer");
            }

            OEAtomBondPath atbPath = dfs.next();
            OEAtomBase at2 = atbPath.getLastNode();
            int j = at2.GetIdx();
            if( j >= i ) continue;

            OEBondBase[] bonds = atbPath.getVertexSequence();
            int mDist = 0;
            for( OEBondBase bd : bonds)
               mDist += multiBondOrdersTimes12[bd.GetIdx()];

            if( multiDMat12[i][j] > mDist )
               multiDMat12[i][j] = mDist;

            if( multiLongDMat12[i][j] < mDist )
               multiLongDMat12[i][j] = mDist;
         }
         dfs.close();
      }
      atIt.delete();

      for( int i=0; i<maxAtIdx+1; i++)
      {  multiDMat12[i][i] = 0;
         multiLongDMat12[i][i] = 0;
         for(int j=i+1; j<maxAtIdx+1; j++)
         {  multiDMat12[i][j] = multiDMat12[j][i];
            multiLongDMat12[i][j] = multiLongDMat12[j][i];
         }
      }


      return multiDMat12;
   }


   /** Multi Bondorder times 12.
    * 8 for aromatic bonds.
    * 12/Bondorder for others
    * This allows for integer values
    */
   private static int getMultiBondOrderTimes12(OEBondBase bd)
   {  if (bd.IsAromatic())
         return 8;
      return 12/bd.GetOrder();
   }


   /** Eq. 8 in Ref 3 */
   private double getBalabanDY(int atIdxi, double atDegree)
   {  OEAtomBase at = ats[atIdxi];
      int z = at.GetAtomicNum();
      int g = Atom.MAIN_GROUP[z];

      double d = atDegree * (1.1191 + 0.016 * z - 0.0537 * g);

      return d;
   }

   /** Eq. Xi form Ref 5, page 23
    * @param multiLongDegrees2 */
   private double getBalabanDX(int atIdx, double atDegree)
   {  OEAtomBase at = ats[atIdx];
      int z = at.GetAtomicNum();
      int g = Atom.MAIN_GROUP[z];

      double d = atDegree * (0.4196 - 0.0078 * z + 0.1567 * g);

      return d;
   }


   public static void main(String ... args)
   {  // run TopologicalIndexTest
   }
}
