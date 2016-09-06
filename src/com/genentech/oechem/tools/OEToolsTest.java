/*
   Copyright 2008-2014 Genentech Inc.

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

import openeye.oechem.OEGraphMol;

import org.testng.annotations.Test;

public class OEToolsTest
{
   @Test()
   public static void TestNeutralize()
   {  OEGraphMol mol = new OEGraphMol();

      OETools.smiToMol(mol,"CC(=O)[O-]");
      OETools.neutralize(mol);
      String s = OETools.molToCanSmi(mol, true);
      assert "CC(=O)O".equals(s) : "neutralize CC(=O)[O-]: " + s;
         
      OETools.smiToMol(mol,"CC[N+H3]");
      OETools.neutralize(mol);
      s = OETools.molToCanSmi(mol, true);
      assert "CCN".equals(s) : "neutralize CC[N+H3]: " + s;
         
      OETools.smiToMol(mol,"CC[N+](C)(C)[H]");
      OETools.neutralize(mol);
      s = OETools.molToCanSmi(mol, true);
      assert "CCN(C)C".equals(s) : "neutralize CC[N+](C)(C)[H]: " + s;
         
      OETools.smiToMol(mol,"CC[N+](C)(C)C");
      OETools.neutralize(mol);
      s = OETools.molToCanSmi(mol, true);
      assert "CC[N+](C)(C)C".equals(s) : "neutralize CC[N+](C)(C)C: " + s;
         
      OETools.smiToMol(mol,"CC[N+](=O)[O-]");
      OETools.neutralize(mol);
      s = OETools.molToCanSmi(mol, true);
      assert "CC[N+](=O)[O-]".equals(s) : "neutralize CC[N+](=O)[O-]: " + s;      
      
      mol.delete();
   }
}
