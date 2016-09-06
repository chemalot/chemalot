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

public class IAAPathGeneratorCharTest
{  String[] smis = new String[] { "C", "C(=O)F", "C1ON1" };
   int[][][] expectedPaths =                                        // int's are ordered ascending, strings are not ordered
            {{ { } },                                               //  { "" }
             { { 1160, 2270 }, { 2260, 30692 }, { 1145, 33761 } },  //  { "-F", "=O"} , {"=C", "=C-F"} , { "-C", "-C=O"}
             { { 752,  1150, 1155, 3826,  38221, 43791},            //  { "-N-O-", "-N-O", "-N", "-O-N-", "-O-N", "-O" }
               { 1145, 1150, 1596, 4670,  32641, 38211},            //  { "-C-N-", "-C-N", "-C", "-N-C-", "-N-C", "-N" }
               { 1145, 1155, 5785, 32646, 43786, 65173}}};          //  { "-O-C-", "-O-C", "-O", "-C-O-", "-C-O", "-C" }

   @Test
   public void testAPG()
   {  OEMolBase mol = new OEGraphMol();
      for( int i=0; i<smis.length; i++)
     {  String smi = smis[i];
        int[][] expectedPath = expectedPaths[i];
        mol.Clear();
        oechem.OEParseSmiles(mol, smi);
        IAAPathGeneratorChar aapg = new IAAPathGeneratorChar(mol, 7);

        int atIdx = 0;
        OEAtomBaseIter atit = mol.GetAtoms();
        while(atit.hasNext())
        {  char[] path = aapg.getAtomPaths(atit.next());

           assert path.length == expectedPath[atIdx].length:
              String.format("Did not get the expected path length %d != %d for %s\n",
                       path.length, expectedPath.length, smi);

           for(int j=0; j<expectedPath[atIdx].length; j++)
              assert path[j] == expectedPath[atIdx][j] :
                 String.format("Path missmatch for %s: %d != %s\n",
                          smi, (int)path[j], expectedPath[atIdx][j]);
           atIdx++;
        }
        atit.delete();
      }
      mol.delete();
   }

}
