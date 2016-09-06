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

import com.aestel.chemistry.openEye.SimComparatorFactory;

public class FPComparatorFact implements SimComparatorFactory<OEMolBase, FPComparator, FPComparator>
{  private final String fpTag;
   private final boolean doMaxTanimoto;

   public FPComparatorFact(boolean doMaxTanimoto, String fpTag)
   {  this.fpTag = fpTag;
      this.doMaxTanimoto = doMaxTanimoto;
   }

   public FPComparator createComparable(OEMolBase baseObject)
   {  if( doMaxTanimoto )
         return new FPMTaniComparator(baseObject, fpTag);
         
      return new FPComparator(baseObject, fpTag);
   }


   public FPComparator createComparator(FPComparator comparator)
   {  return comparator.clone(); //for fingerprints the comparator == the comparable
   }


   public void close()
   {  // nothing to do
   }

}
