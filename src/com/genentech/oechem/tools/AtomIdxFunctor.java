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
package com.genentech.oechem.tools;

import openeye.oechem.OEAtomBase;
import openeye.oechem.OEUnaryAtomBoolFunc;
import openeye.oechem.OEUnaryAtomPred;

/**
 * An OpenEye functor which can be used eg. in for loops to find bonds of a 
 * specific order.
 * 
 * @author albertgo
 *
 */
public class AtomIdxFunctor extends OEUnaryAtomPred {

   private final int requiredIdx;
   
   public AtomIdxFunctor(int idx) {
      requiredIdx = idx;
   }
   

   @Override
   public boolean constCall(OEAtomBase at) {
      return ( at.GetIdx() == requiredIdx ) ;
   }
   
   @Override
   public OEUnaryAtomBoolFunc CreateCopy( ) {
      OEUnaryAtomBoolFunc copy = new AtomIdxFunctor(requiredIdx);
      copy.swigReleaseOwnership();
      return copy;
   }
}
