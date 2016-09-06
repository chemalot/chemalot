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
package com.genentech.chemistry.openEye;

import openeye.oechem.*;

import org.testng.annotations.Test;

public class AAPathGeneratorTest
{  String[] smis = new String[] { "C", "C(=O)F", "C1ON1" };
   String[][][] expectedPaths =
            {{ { "C" } },
             { {"C=O", "C-F", "C"} , {"O=C-F", "O=C", "O"} , {"F-C=O", "F-C", "F"}},
             { { "C-N-O-", "C-N-O", "C-N", "C-O-N-", "C-O-N", "C-O", "C" },
               { "O-C-N-", "O-C-N", "O-C", "O-N-C-", "O-N-C", "O-N", "O" },
               { "N-O-C-", "N-O-C", "N-O", "N-C-O-", "N-C-O", "N-C", "N" }}};

   @Test
   public void testAPG()
   {  OEMolBase mol = new OEGraphMol();
      for( int i=0; i<smis.length; i++)
     {  String smi = smis[i];
        String[][] expectedPath = expectedPaths[i];
        mol.Clear();
        oechem.OEParseSmiles(mol, smi);
        AtomPathGenerator aapg = new AtomPathGenerator(mol, 7);

        int atIdx = 0;
        OEAtomBaseIter atit = mol.GetAtoms();
        while(atit.hasNext())
        {  String[] path = aapg.getAtomPaths(atit.next());
           assert path.length == expectedPath[atIdx].length:
              String.format("Did not get the expected path length %d != %d for %s\n",
                       path.length, expectedPath.length, smi);

           for(int j=0; j<expectedPath[atIdx].length; j++)
              assert path[j].equals(expectedPath[atIdx][j]) :
                 String.format("Path missmatch for %s: %s != %s\n",
                          smi, path[j], expectedPath[atIdx][j]);
           atIdx++;
        }
        atit.delete();
      }
      mol.delete();
   }

}
