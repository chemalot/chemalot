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
package com.genentech.struchk.oeStruchk;

import openeye.oechem.OEGraphMol;
import openeye.oechem.oechem;

import com.aestel.utility.DataFormat;
import com.genentech.oechem.tools.OETools;

/**
 * Just a container for salt data.
 */
class Salt {
   final String code;
   final String name;
   final String canISmiles;
   final String mw;
   final String mf;

   Salt(String smi, String code, String name) {
      this.code = code;
      this.name = name;

      if(smi.length() == 0) {
         // allow for parent salt code
         this.canISmiles = "";
         mw = "0";
         mf = "";
         return;
      }
      
      OEGraphMol mol = new OEGraphMol();
     
      OETools.smiToMol(mol, smi);
      if(mol.NumAtoms() == 0)
         throw new Error(String.format("Invalid smiles: %s code=%s.", smi, code));
      
      canISmiles = OETools.molToCanSmi(mol, true);
      mw = DataFormat.formatNumber(oechem.OECalculateMolecularWeight(mol),"r2");
      mf = oechem.OEMolecularFormula(mol);
      mol.delete();
   }
}
