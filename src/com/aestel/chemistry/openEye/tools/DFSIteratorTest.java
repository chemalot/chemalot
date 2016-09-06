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

package com.aestel.chemistry.openEye.tools;
import openeye.oechem.OEAtomBaseIter;
import openeye.oechem.OEGraphMol;
import openeye.oechem.oechem;

import org.testng.annotations.Test;

import com.genentech.oechem.tools.OETools;

public class DFSIteratorTest
{  @Test
   public void testLinear()
   {  assertPath("C", "C");
      assertPath("FC(Cl)[H]","Cl-C-F C-F F C-F C-Cl C Cl-C-F C-Cl Cl");
      assertPath("c1[nH]ccc1CN1[SiH](F)=[P](Br)1", 
            "c:c:c:c:n: c:c:c:c:n c:c:c:c c:c:c c:c-C-N-P=Si- F-Si=P-N-C-c:c "
           +"Si=P-N-C-c:c Br-P-N-C-c:c P-N-C-c:c F-Si-N-C-c:c Br-P=Si-N-C-c:c "
           +"c:c-C-N-P=Si- P=Si-N-C-c:c Si-N-C-c:c N-C-c:c C-c:c c:c c:c:c:c:n: "
           +"P-N-C-c:c:c:n:c Si-N-C-c:c:c:n:c N-C-c:c:c:n:c C-c:c:c:n:c c:c:c:n:c "
           +"c:c:n:c c:n:c c:n c c:c:c:c:n: c:c:c:c:n c:c:c:n n:c:c-C-N-P=Si- "
           +"F-Si=P-N-C-c:c:n Si=P-N-C-c:c:n Br-P-N-C-c:c:n P-N-C-c:c:n "
           +"F-Si-N-C-c:c:n Br-P=Si-N-C-c:c:n n:c:c-C-N-P=Si- P=Si-N-C-c:c:n "
           +"Si-N-C-c:c:n N-C-c:c:n C-c:c:n c:c:n c:n c:c:c:c:n: c:c:c:c:n "
           +"Si=P-N-C-c:c:c:n Br-P-N-C-c:c:c:n P-N-C-c:c:c:n F-Si-N-C-c:c:c:n "
           +"P=Si-N-C-c:c:c:n Si-N-C-c:c:c:n N-C-c:c:c:n C-c:c:c:n c:c:c:n c:c:n "
           +"c:n n c:c:c:c:n: c:c:c:n:c Si=P-N-C-c:c:n:c Br-P-N-C-c:c:n:c "
           +"P-N-C-c:c:n:c F-Si-N-C-c:c:n:c P=Si-N-C-c:c:n:c Si-N-C-c:c:n:c "
           +"N-C-c:c:n:c C-c:c:n:c c:c:n:c c:n:c c:n c:c:c:c:n: c:c:c:c:n c:c:c:c "
           +"c:c:c-C-N-P=Si- F-Si=P-N-C-c:c:c Si=P-N-C-c:c:c Br-P-N-C-c:c:c "
           +"P-N-C-c:c:c F-Si-N-C-c:c:c Br-P=Si-N-C-c:c:c c:c:c-C-N-P=Si- "
           +"P=Si-N-C-c:c:c Si-N-C-c:c:c N-C-c:c:c C-c:c:c c:c:c c:c c c:c:c:c:n: "
           +"P-N-C-c:c:n:c:c Si-N-C-c:c:n:c:c N-C-c:c:n:c:c C-c:c:n:c:c c:c:n:c:c "
           +"c:c:n:c c:c:n c:c c:c:c:c:n: c:c:c:n:c c:c:c:n c:c:c c:c-C-N-P=Si- "
           +"F-Si=P-N-C-c:c Si=P-N-C-c:c Br-P-N-C-c:c P-N-C-c:c F-Si-N-C-c:c "
           +"Br-P=Si-N-C-c:c c:c-C-N-P=Si- P=Si-N-C-c:c Si-N-C-c:c N-C-c:c C-c:c "
           +"c:c c c:c:c:c:n: c:c:c:n:c c:c:c:n c:c:c c:c c:c:c:c:n: c:c:n:c:c "
           +"c:c:n:c c:c:n c:c c-C-N-P=Si- F-Si=P-N-C-c Si=P-N-C-c Br-P-N-C-c "
           +"P-N-C-c F-Si-N-C-c Br-P=Si-N-C-c c-C-N-P=Si- P=Si-N-C-c Si-N-C-c "
           +"N-C-c C-c c C-c:c:c:n:c: C-c:c:c:n:c C-c:c:c:n C-c:c:c C-c:c C-c:c:c:n:c: "
           +"C-c:c:n:c:c C-c:c:n:c C-c:c:n C-c:c C-c C-N-P=Si- C-N-P=Si-F C-N-P=Si "
           +"Br-P-N-C C-N-P C-N-Si-F Br-P=Si-N-C C-N-P=Si- C-N-Si=P C-N-Si C-N C "
           +"N-C-c:c:c:n:c: N-C-c:c:c:n:c N-C-c:c:c:n N-C-c:c:c N-C-c:c N-C-c:c:c:n:c: "
           +"N-C-c:c:n:c:c N-C-c:c:n:c N-C-c:c:n N-C-c:c N-C-c C-N N-P=Si- "
           +"F-Si=P-N N-P=Si Br-P-N N-P F-Si-N Br-P=Si-N N-P=Si- N-Si=P N-Si N "
           +"Si-N-C-c:c:c:n:c Si-N-C-c:c:c:n Si-N-C-c:c:c Si-N-C-c:c Si-N-C-c:c:n:c:c "
           +"Si-N-C-c:c:n:c Si-N-C-c:c:n Si-N-C-c:c Si-N-C-c C-N-Si N-P=Si- "
           +"Br-P-N-Si P-N-Si N-Si F-Si Br-P=Si Si=P-N-C-c:c:c:n Si=P-N-C-c:c:c "
           +"Si=P-N-C-c:c Si=P-N-C-c:c:n:c Si=P-N-C-c:c:n Si=P-N-C-c:c Si=P-N-C-c "
           +"C-N-P=Si N-P=Si- N-P=Si P=Si Si F-Si-N-C-c:c:c:n F-Si-N-C-c:c:c "
           +"F-Si-N-C-c:c F-Si-N-C-c:c:n:c F-Si-N-C-c:c:n F-Si-N-C-c:c F-Si-N-C-c "
           +"C-N-Si-F F-Si-N-P= Br-P-N-Si-F F-Si-N-P F-Si-N Br-P=Si-F F-Si=P-N-C-c:c:c "
           +"F-Si=P-N-C-c:c F-Si=P-N-C-c:c:n F-Si=P-N-C-c:c F-Si=P-N-C-c C-N-P=Si-F "
           +"F-Si-N-P= F-Si=P-N F-Si=P F-Si F P=Si-N-C-c:c:c:n P=Si-N-C-c:c:c "
           +"P=Si-N-C-c:c P=Si-N-C-c:c:n:c P=Si-N-C-c:c:n P=Si-N-C-c:c P=Si-N-C-c "
           +"C-N-Si=P N-P=Si- N-Si=P F-Si=P P=Si Br-P P-N-C-c:c:c:n:c P-N-C-c:c:c:n "
           +"P-N-C-c:c:c P-N-C-c:c P-N-C-c:c:n:c:c P-N-C-c:c:n:c P-N-C-c:c:n P-N-C-c:c "
           +"P-N-C-c C-N-P F-Si-N-P N-P=Si- P-N-Si N-P P Br-P=Si-N-C-c:c:c "
           +"Br-P=Si-N-C-c:c Br-P=Si-N-C-c:c:n Br-P=Si-N-C-c:c Br-P=Si-N-C-c "
           +"Br-P=Si-N-C Br-P-N-Si= Br-P=Si-N Br-P=Si-F Br-P=Si Br-P-N-C-c:c:c:n "
           +"Br-P-N-C-c:c:c Br-P-N-C-c:c Br-P-N-C-c:c:n:c Br-P-N-C-c:c:n Br-P-N-C-c:c "
           +"Br-P-N-C-c Br-P-N-C Br-P-N-Si-F Br-P-N-Si= Br-P-N-Si Br-P-N Br-P Br");

    }  
   
   @Test
   public void testRing()
   {  assertPath("C#1-O=N#1", "C#N=O- C#N=O C#N C#N=O- C-O=N C-O C "
                             +"C#N=O- N#C-O C-O C#N=O- C#N=O N=O O "
                             +"C#N=O- C-O=N N=O C#N=O- N#C-O C#N N");
      assertPath("O-1=N#C-1", "C#N=O- N#C-O C-O C#N=O- C#N=O N=O O "
                             +"C#N=O- C-O=N N=O C#N=O- N#C-O C#N N "
                             +"C#N=O- C#N=O C#N C#N=O- C-O=N C-O C");
      assertPath("N=1#C-O=1", "C#N=O- C-O=N N=O C#N=O- N#C-O C#N N "
                             +"C#N=O- C#N=O C#N C#N=O- C-O=N C-O C "
                             +"C#N=O- N#C-O C-O C#N=O- C#N=O N=O O");
      assertPath("C-1#N=O-1", "C#N=O- C-O=N C-O C#N=O- C#N=O C#N C "
                             +"C#N=O- N#C-O C#N C#N=O- C-O=N N=O N "
                             +"C#N=O- C#N=O N=O C#N=O- N#C-O C-O O");
      
      assertPath("ClC#1-O=N#1",
           "Cl-C#N=O- Cl-C#N=O Cl-C#N Cl-C#N=O- Cl-C-O=N Cl-C-O C-Cl Cl "
          +"C-Cl C#N=O- C#N=O C#N C#N=O- C-O=N C-O C "
          +"Cl-C-O C#N=O- N#C-O C-O Cl-C#N=O C#N=O- C#N=O N=O O "
          +"Cl-C-O=N C#N=O- C-O=N N=O Cl-C#N C#N=O- N#C-O C#N N");

      assertPath("c1[nH]ccc1", "c:c:c:c:n: c:c:c:c:n c:c:c:c c:c:c c:c "
                              +"c:c:c:c:n: c:c:c:n:c c:c:n:c c:n:c c:n c "
                              +"c:c:c:c:n: c:c:c:c:n c:c:c:n c:c:n c:n "
                              +"c:c:c:c:n: c:c:c:c:n c:c:c:n c:c:n c:n n "
                              +"c:c:c:c:n: c:c:c:n:c c:c:n:c c:n:c c:n "
                              +"c:c:c:c:n: c:c:c:c:n c:c:c:c c:c:c c:c c "
                              +"c:c:c:c:n: c:c:n:c:c c:c:n:c c:c:n c:c "
                              +"c:c:c:c:n: c:c:c:n:c c:c:c:n c:c:c c:c c "
                              +"c:c:c:c:n: c:c:c:n:c c:c:c:n c:c:c c:c "
                              +"c:c:c:c:n: c:c:n:c:c c:c:n:c c:c:n c:c c");
   }
   
   
   static void assertPath(String smi, String pathStr)
   {  String outPath = getPathStr(smi);
      assert pathStr.equals(outPath) : "got: " + outPath + "\texpected: " + pathStr;
   }
   
   static String getPathStr(String smi)
   {  OEGraphMol mol = new OEGraphMol();
      OETools.smiToMol(mol, smi);
      oechem.OESuppressHydrogens(mol);
      StringBuilder sb = new StringBuilder(smi.length()*5);
      
      OEAtomBaseIter aIt = mol.GetAtoms();
      while(aIt.hasNext())
      {  DFSIterator dfs = new DFSIterator(aIt.next(), 7);
         while(dfs.hasNext())
         {  sb.append(DFSIterator.getCannonicalPath(dfs.next())).append(" ");
         }
         dfs.close();
      }
      aIt.delete();
      mol.delete();
      return sb.toString().trim();
   }
}
