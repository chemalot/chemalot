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
 * provides atom types corresponding to th atoms and bonds in the smiles string.
 *
 * @author albertgo
 *
 */
public class SmilesTyper implements AtomTyperInterface, BondTyperInterface
{  /** this is thread safe so only one instance is needed */
   public static final SmilesTyper INSTANCE = new SmilesTyper();

   private static final String[] BONDTypes = {"~","-", "=", "#", "%"};
   private static final String   AROMATICBond = ":";
   private static final String[] ATOMTypes;
   private static final String[] AROMATICAtomTypes;

   static
   {  ATOMTypes = new String[108];
      AROMATICAtomTypes= new String[108];
      ATOMTypes[0] = "*";
      AROMATICAtomTypes[0] = "*";

      for(int i=1; i<108; i++)
      {  String sym = oechem.OEGetAtomicSymbol(i);
         ATOMTypes[i] = sym;
         AROMATICAtomTypes[i] = sym.toLowerCase();
      }
   }
   @Override
   public int getIntType(OEAtomBase at)
   {  int atNum = at.GetAtomicNum();
      if(at.IsAromatic()) atNum += 128;
      return atNum;
   }

   @Override
   public String getType(OEAtomBase at)
   {  int aNum = at.GetAtomicNum();

      if(at.IsAromatic())
         return AROMATICAtomTypes[aNum];
      return ATOMTypes[aNum];
   }

   @Override
   public void initType(OEMolBase mol)
   {  // nothing to do here
   }

   @Override
   public int getIntType(OEBondBase bd)
   {  if(bd.IsAromatic()) return 0;
      return bd.GetOrder();
   }

   @Override
   public String getType(OEBondBase bd)
   {  if(bd.IsAromatic()) return AROMATICBond;
      return BONDTypes[bd.GetOrder()];
   }
}
