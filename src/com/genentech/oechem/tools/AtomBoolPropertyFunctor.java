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
 * An OpenEye functor which can be used to find atoms with a given value for
 * a given property using the {@link OEAtomBase#GetBoolData(int)} method.
 * 
 * @author albertgo
 *
 */
public class AtomBoolPropertyFunctor extends OEUnaryAtomPred {

   private final int requiredTag;
   private final boolean requiredValue;
   
   public AtomBoolPropertyFunctor(int tag, boolean value) {
      requiredTag = tag;
      requiredValue = value;
   }
   
   
   @Override
   public boolean constCall(OEAtomBase at) {
      return ( at.GetBoolData(requiredTag) == requiredValue ) ;
   }
   
   
   @Override
   public OEUnaryAtomBoolFunc CreateCopy( ) {
      OEUnaryAtomBoolFunc copy = new AtomBoolPropertyFunctor(requiredTag, requiredValue);
      copy.swigReleaseOwnership();
      return copy;
   }
}
