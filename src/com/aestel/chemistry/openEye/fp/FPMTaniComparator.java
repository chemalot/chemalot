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

import com.aestel.chemistry.openEye.SimComparator;

import openeye.oechem.OEMolBase;


/**
 * Computes maximumTanimoto = common / (Max(nbits1,nBits2)*2 - common)
 */
public class FPMTaniComparator extends FPComparator
{  public FPMTaniComparator(OEMolBase thisMol, String fpTag)
   {  super(thisMol, fpTag);
   }

   public double similarity(FPMTaniComparator otherFP)
   {   return fp.mtanimoto(otherFP.fp);
   }

   public double similarity(SimComparator<FPComparator> other)
   {  return similarity((FPMTaniComparator) other);
   }
}
