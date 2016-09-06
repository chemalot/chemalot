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

package com.aestel.chemistry.openEye.fp;

import java.util.NoSuchElementException;

import openeye.oechem.*;

import com.aestel.chemistry.openEye.tools.DFSIterator;
import com.aestel.chemistry.openEye.tools.OEAtomBondPath;

/**
 * Generate linear fingerprint codes by dept first search from each atom in the molecule.
 * @author albertgo
 *
 */
public class LinearCodeNameGenerator implements StructureCodeNameIterator
{  private final String type;

   // for rings we are returning two code names; one with stars and one with atom types
   private String secondRingCodeName = null;
   private String matchedCodeName;
   private DFSIterator dfsIterator = null;
   private OEAtomBaseIter atIter;
   private final int depth;
   private final int starAtomDepth;

   /**
    *
    * @param type unique name for this CodeNameType.
    * @param depth Maximum number of bonds in linear codes
    * @param starAtomDepth at depth &gt;= starAtomDepth only * will be used
    *        for atom symbols
    */
   public LinearCodeNameGenerator(String type, int depth, int starAtomDepth)
   {  this.depth = depth;
      this.starAtomDepth = starAtomDepth;
      this.type = type;
   }

   @Override
   public void init(OEMolBase mol)
   {  oechem.OESuppressHydrogens(mol);
      oechem.OEAssignAromaticFlags(mol);
      this.atIter = mol.GetAtoms();
   }

   /**
    * Are there any more bits for this molecule?
    */
   @Override
   public boolean hasNext()
   {  if( matchedCodeName != null ) return true; // last one was not fetched
      if( secondRingCodeName != null)
      {  // we found a ring and are flagging this by adding a second bit
         matchedCodeName = secondRingCodeName;
         secondRingCodeName = null;
         return true;
      }

      if(atIter == null) return false; // no more bits

      if( dfsIterator != null && ! dfsIterator.hasNext() )
      {  dfsIterator.close();
         dfsIterator = null;
      }

      // next atom in molecule
      if(dfsIterator == null)
      {  if(!atIter.hasNext())
         {  atIter.delete();
            atIter = null;
            return false;
         }

         OEAtomBase at = atIter.next();
         dfsIterator = new DFSIterator(at,depth);
      }

      OEAtomBondPath path = dfsIterator.next();
      if(path.getNVertices() >= starAtomDepth)
      {  matchedCodeName = DFSIterator.getCannonicalPath(path)
                        .toString(StarAtomTyper.INSTANCE, SmilesTyper.INSTANCE);
         if(isRing(path))
            secondRingCodeName = DFSIterator.getCannonicalPath(path).toString();
      } else
      {  matchedCodeName = DFSIterator.getCannonicalPath(path).toString();
      }

      return true;
   }

   private static boolean isRing(OEAtomBondPath path)
   {  if(path.getNNodes() != path.getNVertices()) return false;
      OEAtomBase[] ats = path.getNodeSequence();
      OEAtomBase firstAt = ats[0];
      OEAtomBase lastAt  = ats[ats.length-1];
      OEBondBase lastBd = path.getVertexSequence()[path.getNVertices()-1];
      if(firstAt.GetIdx() != lastBd.GetNbr(lastAt).GetIdx()) return false;

      return true;
   }

   @Override
   public String next()
   {  if(! hasNext() )
         throw new NoSuchElementException();

      String tmp = matchedCodeName;
      matchedCodeName = null;
      return tmp;
   }

   @Override
   public void close()
   {  if(dfsIterator != null) dfsIterator.close();
      if(atIter != null)      atIter.Decrement();
   }

   /**
    * @throws UnsupportedOperationException;
    */
   @Override
   public void remove()
   {  throw new UnsupportedOperationException();
   }

   @Override
   public String getType()
   {  return type;
   }
}
