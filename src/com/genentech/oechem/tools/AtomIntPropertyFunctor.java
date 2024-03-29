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

import openeye.oechem.*;

/**
 * An OpenEye functor which can be used to find atoms with a given value for
 * a given property using the {@link OEAtomBase#GetIntData(int)} method.
 *
 * @author albertgo
 *
 */
public class AtomIntPropertyFunctor extends OEUnaryAtomPred {

   private final int requiredTag;
   private final int requiredValue;

   public AtomIntPropertyFunctor(int tag, int value) {
      requiredTag = tag;
      requiredValue = value;
   }


   @Override
   public boolean constCall(OEAtomBase at) {
      return ( at.GetIntData(requiredTag) == requiredValue ) ;
   }


   @Override
   public OEUnaryAtomBoolFunc CreateCopy( ) {
      OEUnaryAtomBoolFunc copy = new AtomIntPropertyFunctor(requiredTag, requiredValue);
      copy.swigReleaseOwnership();
      return copy;
   }

   public static OEAtomBase getFirst(OEMolBase mol, int tag, int value)
   {  AtomIntPropertyFunctor aipf = new AtomIntPropertyFunctor(tag, value);
      OEAtomBaseIter atIt = mol.GetAtoms(aipf);
      OEAtomBase res = null;
      if( atIt.hasNext() )
         res = atIt.next();
      atIt.delete();
      aipf.delete();

      return res;
   }
}
