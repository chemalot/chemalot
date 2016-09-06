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

import com.aestel.chemistry.openEye.fp.AtomTyperInterface;
import com.aestel.chemistry.openEye.fp.BondTyperInterface;
import com.aestel.chemistry.openEye.fp.SmilesTyper;

import openeye.oechem.OEAtomBase;
import openeye.oechem.OEBondBase;

public class OEAtomBondPath extends Path<OEAtomBase,OEBondBase>
{  private static final SmilesTyper SMILESTyper = SmilesTyper.INSTANCE;
   
   public OEAtomBondPath(OEAtomBase[] nodes, OEBondBase[] vertices)
   {  super(nodes, vertices);
   }

   
   public String toString()
   {  return toString(SMILESTyper, SMILESTyper);
   }
   
   /**
    * The atTyper and bdTyper have to be initialized with the correct mol.
    */
   public String toString(AtomTyperInterface atTyper, BondTyperInterface bdTyper)
   {  OEAtomBase[] nodes = getNodeSequence();
      OEBondBase[] vertices = getVertexSequence();
      
      StringBuilder sb = new StringBuilder((vertices.length + nodes.length)*2);
      int nIdx = 0;
      for(OEAtomBase node : nodes)
      {  sb.append(atTyper.getType(node));
         if(nIdx < vertices.length)
            sb.append(bdTyper.getType(vertices[nIdx++])); 
      }
      while(nIdx < vertices.length)
         sb.append(bdTyper.getType(vertices[nIdx++]));
      
      return sb.toString();
   }
}
