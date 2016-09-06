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
package com.genentech.chemistry.openEye.cats;

import openeye.oechem.OEAtomBase;
import openeye.oechem.OESubSearch;

public class SSTyper implements AtomTyperInterface
{  private final OESubSearch ATSS;
   private final int idx;
   private final String name;

   public SSTyper(String smarts, int idx, String name)
   {  this.idx =idx;
      this.name = name;
      ATSS = new OESubSearch(smarts);
   }


   @Override
   public boolean isType(OEAtomBase at)
   {  return ATSS.AtomMatch(at);
   }

   @Override
   public int getTypeIdx()
   {  return idx;
   }


   @Override
   public String getTypeName()
   {  return name;
   }


   @Override
   public void close()
   {  ATSS.delete();
   }
}
