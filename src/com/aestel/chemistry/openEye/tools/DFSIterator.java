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

import openeye.oechem.*;

import com.aestel.chemistry.openEye.fp.AtomTyperInterface;
import com.aestel.chemistry.openEye.fp.BondTyperInterface;
import com.aestel.chemistry.openEye.fp.SmilesTyper;
import com.genentech.oechem.tools.OETools;

/**
 * Depth first Search implementation.
 *
 * Each call to next() will return the next lienar path that originates at the
 * center atom passed to the constuctor.
 *
 * @author albertgo
 *
 */
public class DFSIterator implements Iterator<OEAtomBondPath> {
   private final boolean[] atVisited;
   private final boolean[] bdVisited;
   private final Deque<OEAtomBase> atPath;
   private final Deque<OEBondBase> bdPath;
   private final Deque<OEBondBaseIter> atBdIterPath;
   private OEAtomBondPath currentPath;
   private int maxBonds;

   public DFSIterator(OEAtomBase centerAt, int maxBonds)
   {  this.maxBonds = maxBonds;
      atVisited = new boolean[ centerAt.GetParent().GetMaxAtomIdx() ];
      bdVisited = new boolean[ centerAt.GetParent().GetMaxBondIdx() ];
      atPath = new LinkedList<OEAtomBase>();
      bdPath = new LinkedList<OEBondBase>();

      /** list of bond itterators around the atoms alonng the current path*/
      atBdIterPath = new LinkedList<OEBondBaseIter>();
      // LinkedList seems to be about 3% faster than ArrayDeque in java 1.6
//      atPath = new ArrayDeque<OEAtomBase>(Math.min(maxBonds+1,at.GetParent().NumAtoms()));
//      bdPath = new ArrayDeque<OEBondBase>(Math.min(maxBonds,at.GetParent().NumBonds()));
//      atBdIterPath = new ArrayDeque<OEBondBaseIter>(Math.min(maxBonds,at.GetParent().NumBonds()));

      atPath.addLast(centerAt);
      atBdIterPath.addLast(centerAt.GetBonds());

      atVisited[centerAt.GetIdx()] = true;
   }


   @Override
   public boolean hasNext() {
      if( currentPath != null ) return true;

      if(atBdIterPath.size() == 0) return false;


      if( ! nextPath() ) return false;

      OEAtomBase[] atoms = atPath.toArray(new OEAtomBase[atPath.size()]);
      OEBondBase[] bonds = bdPath.toArray(new OEBondBase[bdPath.size()]);

      currentPath = new OEAtomBondPath(atoms, bonds);
      return true;
   }


   @Override
   public OEAtomBondPath next()
   {  if(! hasNext() ) throw new NoSuchElementException();

      OEAtomBondPath tmp = currentPath;
      currentPath = null;
      return tmp;
   }


   private boolean nextPath()
   {  if(atBdIterPath.size() == 0) return false;

      if( atPath.size() > atBdIterPath.size() )
      {  // no iterator for this bond => we are at end
         OEAtomBase at = atPath.pollLast();
         OEBondBase bd = bdPath.pollLast();
         // mark as reusable since we might get here from an other direction of
         // a circle
         atVisited[at.GetIdx()] = false;
         bdVisited[bd.GetIdx()] = false;
      }

      if( bdPath.size() == atPath.size() )  // we found ring in last cycle and have dangling bond
         bdVisited[ bdPath.pollLast().GetIdx() ] = false;

      assert atPath.size() == atBdIterPath.size();
      assert bdPath.size()+1 == atPath.size();

      // take next bond from last atom
      OEBondBaseIter bIt = atBdIterPath.peekLast();
      while(bIt.hasNext())
      {  OEBondBase bd = bIt.next();
         if( bdVisited[bd.GetIdx()] ) continue;

         bdPath.addLast(bd);
         bdVisited[bd.GetIdx()] = true;

         // get atom across new bond from last atom
         OEAtomBase lastAt = atPath.peekLast();
         OEAtomBase nextAt = bd.GetNbr(lastAt);
         if( atVisited[nextAt.GetIdx()] ) return true;

         atPath.addLast(nextAt);
         atVisited[nextAt.GetIdx()] = true;

         if( nextAt.GetExplicitDegree() < 2) return true;
         if( bdPath.size() == maxBonds ) return true;

         atBdIterPath.addLast(nextAt.GetBonds());
         nextPath(); // depth first recurse one more level

         // do not delete bIt! it is owned by the instance and reused.
         return true;
      }
      atBdIterPath.pollLast().delete(); // this is the currently used bIt

      return true;
   }


   public void close()
   {  for(OEBondBaseIter bIt : atBdIterPath)
         bIt.delete();
   }


   /**
    * throws UnsupportedOperationException
    */
   @Override
   public void remove()
   {  throw new UnsupportedOperationException();
   }

   public static OEAtomBondPath getCannonicalPath(OEAtomBondPath inPath)
   {  return getCannonicalPath(inPath, SmilesTyper.INSTANCE, SmilesTyper.INSTANCE);
   }

   /**
    * get lexically smalles representation of this path.
    */
   public static OEAtomBondPath getCannonicalPath(OEAtomBondPath inPath,
               AtomTyperInterface atomTyper, BondTyperInterface bondTyper)
   {  // these are copies not the originals
      OEAtomBase[] atoms = inPath.getNodeSequence();
      OEBondBase[] bonds = inPath.getVertexSequence();
      if( atoms.length == bonds.length )
      {  // this is a cycle ending in a bond
         OEBondBase lastBd = bonds[bonds.length-1];
         OEAtomBase lastAt = atoms[atoms.length-1];

         // find index of first ring atom in path, this is the same as the atom
         // bound to the last bond.
         int firstRingAtOEIdx = lastBd.GetNbr(lastAt).GetIdx();
         int firstRingAtIdx = 0;
         while(atoms[firstRingAtIdx].GetIdx() != firstRingAtOEIdx)
            firstRingAtIdx++;

         // we need to compare the path after reversing the atoms in the ring
         // C-C1=N:O-1 to C-C1-O:N=1

         if(firstRingAtIdx == 0)
         {  //@TODO for pure rings (firstRingAtIdx == 0) we should generate
            // path starting at all atoms in both directions and choose only
            // the lexically smallest one.
            return getCannonicalRing(atomTyper, bondTyper, atoms,bonds);
         }

         reverse(atoms, firstRingAtIdx+1, atoms.length);
         reverse(bonds, firstRingAtIdx, bonds.length);
      } else
      {  reverse(atoms, 0, atoms.length);
         reverse(bonds, 0, bonds.length);
      }

      // compare by lexical comparison of strings
      String inPathStr = inPath.toString(atomTyper, bondTyper);
      OEAtomBondPath newPath = new OEAtomBondPath(atoms,bonds);
      String newPathStr = newPath.toString(atomTyper, bondTyper);
      if( inPathStr.compareTo(newPathStr) > 0 ) return newPath;
      return inPath;
   }


   /**
    * Return the path which gives the lexically smalles path around the ring.
    * C1-N=O#1 -> N1=O#C-1 -> O1#C-N=1   reverse
    *   -> C1#O=N-1 -> O1=N-C#1 -> N1-C#O=1
    * @param atoms  a list of atoms in the sequence of the ring
    * @param bonds a list of bonds in the sequence of the ring
    */
   private static OEAtomBondPath getCannonicalRing(
         AtomTyperInterface atomTyper, BondTyperInterface bondTyper,
         OEAtomBase[] atoms, OEBondBase[] bonds)
   {  String smallestPath = new OEAtomBondPath(atoms, bonds).toString(atomTyper, bondTyper);
      OEAtomBase[] smallPathAts = atoms.clone();
      OEBondBase[] smallPathBds = bonds.clone();
      for(int i=1; i<atoms.length; i++)
      {  rotate(atoms);
         rotate(bonds);
         String nPath = new OEAtomBondPath(atoms, bonds).toString(atomTyper, bondTyper);
//System.err.println(nPath + " " + smallestPath);
         if(nPath.compareTo(smallestPath) < 0)
         {  smallestPath = nPath;
            System.arraycopy(atoms, 0, smallPathAts, 0, atoms.length);
            System.arraycopy(bonds, 0, smallPathBds, 0, bonds.length);
         }
      }
      // reverse ring direction
      reverse(atoms,1, atoms.length);
      reverse(bonds, 0, bonds.length);

      for(int i=0; i<atoms.length; i++)
      {  rotate(atoms);
         rotate(bonds);
         String nPath = new OEAtomBondPath(atoms, bonds).toString(atomTyper, bondTyper);
//System.err.println(nPath + " " + smallestPath);
         if(nPath.compareTo(smallestPath) < 0)
         {  smallestPath = nPath;
            System.arraycopy(atoms, 0, smallPathAts, 0, atoms.length);
            System.arraycopy(bonds, 0, smallPathBds, 0, bonds.length);
         }
      }

      return new OEAtomBondPath(smallPathAts, smallPathBds);
   }

   /*
    * Move elements of the array one element to higher numbers putting the last
    * element into the first position.
    *
    *  @return a reference to the input array, rotation is done inplace.
    */
   private static <T> T[] rotate(T [] arr)
   {  T last = arr[arr.length-1];
      System.arraycopy(arr, 0, arr, 1, arr.length-1);
      arr[0] = last;

      return arr;
   }

   /**
    * returns reference to input after reversing the order of the elements
    * from start (inclusive) to end (exclusive)
    */
   private static <T> T[] reverse(T[] array, int start, int end)
   {  // exchange the first and last
      for (int left=start, right=end-1; left<right; left++, right--)
      {  T temp      = array[left];
         array[left] = array[right];
         array[right] = temp;
      }
      return array;
   }


   public static void main(String...args)
   {  OEGraphMol mol = new OEGraphMol();
      OETools.smiToMol(mol, "c1[nH]ccc1CN1[SiH](F)=[P](Br)1");

      OEAtomBaseIter aIt = mol.GetAtoms();
      while(aIt.hasNext())
      {  DFSIterator dfs = new DFSIterator(aIt.next(), 7);
         while(dfs.hasNext())
         {  OEAtomBondPath atBPath = dfs.next();

            System.out.printf("%s\t %s\n",atBPath, getCannonicalPath(atBPath));
         }
         dfs.close();
      }
      aIt.delete();
   }
}
