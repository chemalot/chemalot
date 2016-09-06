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

import openeye.oechem.OEMolBase;
import openeye.oechem.oechem;

import com.aestel.chemistry.openEye.SimComparator;

public class FPComparator implements SimComparator<FPComparator>, Cloneable
{  protected final LongFingerprint fp;

   public FPComparator(OEMolBase thisMol, String fpTag)
   {  String fpStr = oechem.OEGetSDData(thisMol, fpTag);
      this.fp = new LongFingerprint(fpStr);
   }

   @Override
   public double similarity(FPComparator otherFP)
   {   return fp.tanimoto(otherFP.fp);
   }

   @Override
   public double similarity(SimComparator<FPComparator> other)
   {  return similarity((FPComparator) other);
   }

   /**
    * Default implementation to be improved.
    */
   @Override
   public double similarity(SimComparator<FPComparator> other, double minSim)
   { return similarity(other);
   }

   @Override
   public void close()
   {  // nothing to do
   }

   @Override
   public FPComparator clone()
   {  FPComparator newFP;
      try
      {  newFP = (FPComparator) super.clone();
      } catch (CloneNotSupportedException e)
      {   throw new Error(e);
      }
      return newFP;
   }
}
