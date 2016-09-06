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

import java.util.Iterator;

import com.aestel.chemistry.openEye.fp.tools.GeneralHashFunctionLibrary;

/**
 * Use a Hash function to hash structureCodenames to their position in a bitString.
 * If addNewCodes is false this is not mutable.
 *
 * @author albertgo
 *
 */
public class HashStructureCodeMapper implements StructureCodeMapper
{  private final int numIdexes;
   private final int minIdx;


   /**
    * A {@link StructureCodeMapper} which derives the position from a hash function
    *
    */
   public HashStructureCodeMapper(int minIdx, int numIdexes)
   {  this.numIdexes = numIdexes;
      this.minIdx = minIdx;
   }

   /**
    * returns hashValue between minIdx (incl.) and minIdx+numIndexes (excl.)
    * for codeName.
    *
    * Currently type is ignored.
    */
   @Override
   public int getIndex(String type, String codeName)
   {  return GeneralHashFunctionLibrary.toInt(
         GeneralHashFunctionLibrary.BKDRHash(codeName), numIdexes)+minIdx;
   }

   /**
    * Returns the highest position number to be returned by {@link #getIndex}.
    */
   @Override
   public int getMaxIdx()
   {  return minIdx+numIdexes-1;
   }

   @Override
   public void close()
   {  // nothing to do here
   }

   @Override
   public int getMinIdx()
   {  return minIdx;
   }

   @Override
   public String getName(int idx)
   {  return Integer.toString(idx);
   }

   @Override
   public Iterator<String[]> getIterator()
   {  // TODO implement iterator from minIdx to maxIdx
      return null;
   }
}
