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
package com.genentech.chemistry.openEye.cats;

import java.util.ArrayList;
import java.util.Arrays;

import openeye.oechem.OEAtomBase;
import openeye.oechem.OEMatchBase;
import openeye.oechem.OEMatchPairAtomIter;

public class Path
{  private final int[] atomIdx;
   private final int hashCode;

   public Path(OEMatchBase m)
   {  ArrayList<Integer> atList = new ArrayList<Integer>();
      OEMatchPairAtomIter atIt = m.GetAtoms();
      while(atIt.hasNext())
      {  int at = atIt.next().getTarget().GetIdx();
         atList.add(at);
      }
      atIt.delete();

      atomIdx = new int[atList.size()];
      for(int i=0; i<atList.size(); i++)
         atomIdx[i] = atList.get(i);

      Arrays.sort(atomIdx);

      hashCode = Arrays.hashCode(atomIdx);
   }

   public Path(OEAtomBase[] ats)
   {  atomIdx = new int[ats.length];
      for(int i=0; i<ats.length; i++)
         atomIdx[i] = ats[i].GetIdx();

      Arrays.sort(atomIdx);

      hashCode = Arrays.hashCode(atomIdx);
   }

   @Override
   public int hashCode()
   {  return hashCode;
   }

   @Override
   public boolean equals(Object other)
   {  if( other == this )  return true;
      if( ! (other instanceof Path)) return false;

      Path p2 = (Path) other;
      return Arrays.equals(this.atomIdx, p2.atomIdx);
   }
}
