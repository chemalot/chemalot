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

import com.aestel.chemistry.openEye.fp.SmilesTyper;

/**
 * Generate linear path starting from the atoms in the mol passed to the constructor.
 *
 * @author albertgo
 *
 */
public final class AtomPathGenerator
{  private final int maxBonds;
   private final String[] atomSymbol;
   private final String[] bondSymbol;
   private final StringBuilder currentPath;
   private final List<String> pathList = new ArrayList<String>();
   private int nBondsDepth;

   /**
    * @param maxBonds maximum number of bonds each linear path may contain.
    */
   public AtomPathGenerator(OEMolBase mol, int maxBonds)
   {  this.maxBonds = maxBonds;
      atomSymbol = new String[mol.GetMaxAtomIdx()];
      bondSymbol = new String[mol.GetMaxBondIdx()];

      // store array of atom symbols
      OEAtomBaseIter aIt = mol.GetAtoms();
      while( aIt.hasNext() )
      {  OEAtomBase a = aIt.next();
         atomSymbol[a.GetIdx()] = SmilesTyper.INSTANCE.getType(a);
      }
      aIt.delete();

      // store array of bond symbols
      OEBondBaseIter bIt = mol.GetBonds();
      while(bIt.hasNext())
      {  OEBondBase b = bIt.next();
         bondSymbol[b.GetIdx()] = SmilesTyper.INSTANCE.getType(b);
      }
      bIt.delete();

      currentPath = new StringBuilder(maxBonds * 3 + 2);
   }


   /**
    * Get all path starting from at.
    */
   public String[] getAtomPaths(OEAtomBase at)
   {  currentPath.setLength(0);
      pathList.clear();

      boolean[] atomVisited = new boolean[atomSymbol.length];
      boolean[] bondVisited = new boolean[bondSymbol.length];
      Arrays.fill(atomVisited, false);
      Arrays.fill(bondVisited, false);

      nBondsDepth = 0;

      // start recursive Depth First Search
      addPath(at, atomVisited, bondVisited);

      //System.err.println(pathList.size());
      return pathList.toArray(new String[pathList.size()]);
   }

   /**
    * Recursively search for paths starting from at.
    *
    * Fill paths into pathList;
    * @param bondVisited boolean array marking atoms visited so far on this path
    * @param atomVisited boolean array marking bonds visited so far on this path
    *
    * @return String of path which includes at least this atom (at).
    */
   private String addPath(OEAtomBase at, boolean[] atomVisited, boolean[] bondVisited)
   {  int aIdx = at.GetIdx();

      if( atomVisited[aIdx] )
      {  String path = currentPath.toString();
         pathList.add(path);
         return path;
      }
      atomVisited[aIdx] = true;

      currentPath.append(atomSymbol[aIdx]);
      String longerPath = null;
      int len = currentPath.length();

      if( nBondsDepth < maxBonds )
      {  nBondsDepth++;
         OEBondBaseIter bIt = at.GetBonds();
         while( bIt.hasNext() )
         {  OEBondBase b = bIt.next();
            int bIdx = b.GetIdx();
//System.err.printf("bidx=%d neigh=%s\n", bIdx, atomSymbol[b.GetNbr(at).GetIdx()]);
            // avoid going back bonds we have traced before
            if( bondVisited[bIdx] ) continue;

            boolean[] nextBondVisited = bondVisited.clone();
            nextBondVisited[bIdx] = true;

            currentPath.append(bondSymbol[bIdx]);

            // keeping longerPath allows java to reuse memory in substrings
            longerPath = addPath(b.GetNbr(at), atomVisited.clone(), nextBondVisited);
            currentPath.setLength(len);
         }
         bIt.delete();
         nBondsDepth--;
      }

      if( longerPath == null ) longerPath = currentPath.toString();
      String path = longerPath.substring(0, len);
      pathList.add(path);

      return longerPath;
   }
}
