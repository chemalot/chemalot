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

import openeye.oechem.*;

/**
 * This AtomTyper types all atoms as equal "*"
 * @author albertgo
 *
 */
public class StarAtomTyper implements AtomTyperInterface
{  /** this is thread safe so only one instance is needed */
   public static final StarAtomTyper INSTANCE = new StarAtomTyper();

   @Override
   public int getIntType(OEAtomBase at)
   {  return 0;
   }

   @Override
   public String getType(OEAtomBase at)
   {  return "*";
   }

   @Override
   public void initType(OEMolBase mol)
   {  // nothing to do here
   }
}
