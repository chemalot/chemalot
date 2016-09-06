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

import java.util.*;

import openeye.oechem.*;

import com.aestel.math.Prime;

/**
 * Modified algorithm to compute AtomAtomPath similarity.
 *
 * The algorithm to compute the similarity of two atoms is changed as follows:
 *    pathsInCommon / max(nPath(AtomMolA), nPath(AtomMolB))
 *
 * This should prevent the situation where a small substituent results in a higher
 * similarity than a big substituent.
 *
 * This implementation uses a BitSet to store the paths of an atom and an int
 * to encode a single path.
 * @author albertgo
 *
 */
public class IAAPathGeneratorFP
{  private static final int MAXAtomType = Prime.getPrimeLargerEqThan(217);
   private static final int MAXBondType = Prime.getPrimeLargerEqThan(5);
   private static final int MAXAtomNum  = 108;
   private static final int NBITS = 256; // needs to be power of 2
   private final int maxBonds;
   private int[] bondTypes;
   private int[] atomTypes;
   private final List<Integer> pathList = new ArrayList<Integer>();
   private int nBondsDepth;
   private int currentPath;
   private final boolean[] atomVisited;
   private final boolean[] bondVisited;

   IAAPathGeneratorFP(OEMolBase mol, int maxBonds)
   {  this.maxBonds = maxBonds;
      atomTypes = new int[mol.GetMaxAtomIdx()+1];
      bondTypes = new int[mol.GetMaxBondIdx()+1];
      atomVisited = new boolean[atomTypes.length];
      bondVisited = new boolean[bondTypes.length];
      Arrays.fill(atomVisited, false);
      Arrays.fill(bondVisited, false);

      OEAtomBaseIter aIt = mol.GetAtoms();
      while( aIt.hasNext() )
      {  OEAtomBase a = aIt.next();
         int aNum = a.GetAtomicNum();
         if( a.IsAromatic() ) aNum = MAXAtomNum + aNum;
         atomTypes[a.GetIdx()] = aNum;
      }
      aIt.delete();

      OEBondBaseIter bIt = mol.GetBonds();
      while(bIt.hasNext())
      {  OEBondBase b = bIt.next();
         int bNum=b.GetOrder();
         if( b.IsAromatic() ) bNum = 4;
         bondTypes[b.GetIdx()] = bNum;
      }
      bIt.delete();
   }

   /**
    * Get all path starting from at.
    */
   public BitSet getAtomPaths(OEAtomBase at)
   {  pathList.clear();

      // do not include this atom so that we can make Cl more like F later
      // loop over neighbors
      int atIdx = at.GetIdx();
      atomVisited[atIdx] = true;
      nBondsDepth = 1;
      OEBondBaseIter bdIt = at.GetBonds();
      while( bdIt.hasNext() )
      {  OEBondBase bd = bdIt.next();
         int bdIdx = bd.GetIdx();
         currentPath = bondTypes[bdIdx]*MAXAtomType;
         bondVisited[bdIdx] = true;
         OEAtomBase nextAt = bd.GetNbr(at);

         addPath(nextAt);
         bondVisited[bdIdx] = false;
      }
      bdIt.delete();
      atomVisited[atIdx] = false;

      Collections.sort(pathList);

      BitSet fp = new  BitSet(NBITS);
      if( pathList.size() == 0 ) return fp;

      int mask = NBITS-1;
      int lastVal = pathList.get(0)-1;
StringBuilder sb1 = new StringBuilder(200);
StringBuilder sb2 = new StringBuilder(200);
      Random r = new Random();
      boolean newLastVal = true;
      for(int i=0; i<pathList.size(); i++)
      {  int val = pathList.get(i);
         if( val == lastVal )
         {  if( newLastVal  ) r.setSeed(val);
            val = r.nextInt();
            newLastVal = false;
         } else
         {  lastVal = val;
            newLastVal = true;
         }

         int bit = val & mask;
sb1.append(lastVal).append(",\t");
sb2.append(bit).append(",\t\t");
         fp.set(bit);
      }
//System.err.println(sb1);
//System.err.println(sb2);
//System.err.printf("nPath=%d nBits=%d\n",pathList.size(),fp.cardinality());
      return fp;
   }

   /**
    * Recursively search for paths starting from at.
    *
    * Fill paths into pathList;
    * @param bondVisited boolean array marking atoms visited so far on this path
    * @param atomVisited boolean array marking bonds visited so far on this path
    * @param depth depth of current atom from start
    *
    * @return int of path which starts from this atom (at) including it.
    */
   private void addPath(OEAtomBase at)
   {  int aIdx = at.GetIdx();

      atomVisited[aIdx] = true;

      int previousCurrentPath = currentPath;
      // unique hashing: make space in the left of currentPath for next bond
      currentPath = (currentPath + atomTypes[aIdx]) * MAXBondType;
      int thisCurrentPath = currentPath;

      if( nBondsDepth < maxBonds )
      {  nBondsDepth++;
         OEBondBaseIter bIt = at.GetBonds();
         while( bIt.hasNext() )
         {  OEBondBase b = bIt.next();
            int bIdx = b.GetIdx();
//System.err.printf("bidx=%d neigh=%d\n", bIdx, atomTypes[b.GetNbr(at).GetIdx()]);
            // avoid going back bonds we have traced before
            if( bondVisited[bIdx] ) continue;

            // unique hashing: make space in the left of currentPath for next atom
            currentPath = (thisCurrentPath + bondTypes[bIdx]) * MAXAtomType;

            OEAtomBase nextAt = b.GetNbr(at);
            if( atomVisited[nextAt.GetIdx()] )
               pathList.add(currentPath); // ring: path contains bond but not atom
            else
            {  bondVisited[bIdx] = true;
               addPath(nextAt);
               bondVisited[bIdx] = false;
            }
         }
         bIt.delete();
         nBondsDepth--;
      }

      pathList.add(thisCurrentPath);
      currentPath = previousCurrentPath;
      atomVisited[aIdx] = false;
   }

   public int[] getAtomTypes()
   {  return atomTypes;
   }

   public static final int atomTypeToAtomNum(int atomType)
   {  if( atomType > MAXAtomNum ) return atomType - MAXAtomNum;
      return atomType;
   }

   public static boolean isAromatic(int atomType)
   {  return atomType > MAXAtomNum;
   }
}
