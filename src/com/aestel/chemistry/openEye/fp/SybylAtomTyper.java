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

import openeye.oechem.OEAtomBase;
import openeye.oechem.OEMolBase;
import openeye.oechem.oechem;

/**
 * Type atoms with Sybyl atom types.
 * @author albertgo
 *
 */
public class SybylAtomTyper implements AtomTyperInterface
{  @Override
public int getIntType(OEAtomBase at)
   {  return at.GetIntType();
   }

   @Override
   public String getType(OEAtomBase at)
   {  return oechem.OETriposTypeName(at.GetIntType());
   }

   @Override
   public void initType(OEMolBase mol)
   {  oechem.OETriposAtomTypes(mol);
   }

}
