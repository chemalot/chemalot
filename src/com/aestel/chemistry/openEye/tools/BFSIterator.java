/*
   Copyright 2006-2014 Man-Ling Lee & Alberto Gobbi

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

Contact: aestelSW@gmail.com
*/

package com.aestel.chemistry.openEye.tools;
import java.util.*;

import openeye.oechem.OEAtomBase;
import openeye.oechem.OEAtomBaseIter;

/**
 * Performs a breath first search around and center atom and provides an
 * iterator access to the atoms sourounding the center.
 *
 * Each time the {@link #next()} method is call an additional surrounding
 * atom si returned. the current depth is returned by {@link #getDepth()}
 *
 * @author albertgo
 *
 */
public class BFSIterator implements Iterator<OEAtomBase> {
   private final boolean[] visited;
   private final Deque<OEAtomBase> queue;
   private final List<OEAtomBase> sphereAtoms;
   private int depth;
   private int maxDepth;

   public BFSIterator(OEAtomBase centerAt, int maxDepth) {
      visited = new boolean[ centerAt.GetParent().GetMaxAtomIdx() ];
      queue = new ArrayDeque<OEAtomBase>(visited.length);
      sphereAtoms = new ArrayList<OEAtomBase>(visited.length);
      this.maxDepth = maxDepth;

      queue.addLast(centerAt);
      visited[centerAt.GetIdx()] = true;
      depth = 0;
      sphereAtoms.add(centerAt);
   }

   @Override
   public boolean hasNext() {
      if( queue.size() > 0 ) return true;

      if( depth >= maxDepth ) return false;

      for(OEAtomBase lastSphereAt : sphereAtoms ) {
         OEAtomBaseIter aIt = lastSphereAt.GetAtoms();
         while(aIt.hasNext()) {
            OEAtomBase at = aIt.next();
            if(! visited[at.GetIdx()] ) {
               queue.addLast(at);
               visited[at.GetIdx()] = true;
            }
         }
         aIt.delete();
      }
      depth++;
      if( queue.size() == 0 ) return false;

      sphereAtoms.clear();
      for( OEAtomBase at : queue)
         sphereAtoms.add(at);

      return true;
   }

   @Override
   public OEAtomBase next() {
      if(! hasNext() ) throw new NoSuchElementException();

      return queue.pollLast();
   }

   /**
    * throws UnsupportedOperationException
    */
   @Override
   public void remove() {
      throw new UnsupportedOperationException();
   }

   /**
    * Depth of atom processed by the last calls to {@link #hasNext()}
    * or {@link #next()} whichever was last.
    */
   public int getDepth() {
      return depth;
   }

}
