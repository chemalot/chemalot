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

import org.jdom.Element;

import com.aestel.utility.MessageList;
import com.genentech.struchk.oeStruchk.OEStruchk.StructureFlag;

public class HydrogenRemover extends AbstractStructureCheck {

   public HydrogenRemover(Element elem) {
      super(elem);
   }

   @Override
   public boolean checkStructure(OEGraphMol in, StructureFlag inStereo,
                                 MessageList msgs) {
      oechem.OESuppressHydrogens(in,false,false,true);
      oechem.OEPerceiveChiral(in);

      return true;
   }

   @Override
   public void delete() { /* nothing to do */ }
}
