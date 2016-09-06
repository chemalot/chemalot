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
import java.util.List;

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
 * This implementation uses a long to store a path.
 * @author albertgo
 *
 */
public class IAAPathGenerator
{  private static final long MAXAtomType = Prime.getPrimeLargerEqThan(217);
   private static final long MAXBondType = Prime.getPrimeLargerEqThan(5);
   private static final int  MAXAtomNum  = 108;
   private final int maxBonds;
   private long[] bondTypes;
   private long[] atomTypes;
   private final List<Long> pathList = new ArrayList<Long>();
   private int nBondsDepth;
   private long currentPath;
   private final boolean[] atomVisited;
   private final boolean[] bondVisited;

   IAAPathGenerator(OEMolBase mol, int maxBonds)
   {  this.maxBonds = maxBonds;
      atomTypes = new long[mol.GetMaxAtomIdx()+1];
      bondTypes = new long[mol.GetMaxBondIdx()+1];
      atomVisited = new boolean[atomTypes.length];
      bondVisited = new boolean[bondTypes.length];
      Arrays.fill(atomVisited, false);
      Arrays.fill(bondVisited, false);

      OEAtomBaseIter aIt = mol.GetAtoms();
      while( aIt.hasNext() )
      {  OEAtomBase a = aIt.next();
         long aNum = a.GetAtomicNum();
         if( a.IsAromatic() ) aNum = MAXAtomNum + aNum;
         atomTypes[a.GetIdx()] = aNum;
      }
      aIt.delete();

      OEBondBaseIter bIt = mol.GetBonds();
      while(bIt.hasNext())
      {  OEBondBase b = bIt.next();
         long bNum=b.GetOrder();
         if( b.IsAromatic() ) bNum = 4;
         bondTypes[b.GetIdx()] = bNum;
      }
      bIt.delete();
   }

   /**
    * Get all path starting from at.
    */
   public long[] getAtomPaths(OEAtomBase at)
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

      long[] paths = new long[pathList.size()];
      for( int i=0; i<pathList.size(); i++ )
         paths[i] = pathList.get(i);

      //System.err.println(paths.length);
      Arrays.sort(paths);
      return paths;
   }

   /**
    * Recursively search for paths starting from at.
    *
    * Fill paths into pathList;
    * @param bondVisited boolean array marking atoms visited so far on this path
    * @param atomVisited boolean array marking bonds visited so far on this path
    * @param depth depth of current atom from start
    *
    * @return long of path which starts from this atom (at) including it.
    */
   private void addPath(OEAtomBase at)
   {  int aIdx = at.GetIdx();

      atomVisited[aIdx] = true;

      long previousCurrentPath = currentPath;
      // unique hashing: make space in the left of currentPath for next bond
      currentPath = (currentPath + atomTypes[aIdx]) * MAXBondType;
      long thisCurrentPath = currentPath;

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

   public long[] getAtomTypes()
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
