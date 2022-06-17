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
package com.genentech.chemistry.tool.align;

import java.util.ArrayList;
import java.util.List;

import com.genentech.chemistry.tool.align.SDFAlign.OUTType;

import openeye.oechem.*;

class McsAlign extends AbstractAlign
{
   private List<OEMCSSearch> mcsss = new ArrayList<OEMCSSearch>();

   public McsAlign(List<OEMol> refmols, OUTType oMolType, String rmsdTag, String atomDevTag,
                   String[] atomInfoTags, boolean doOptimize, boolean doMirror, 
                   int atomExpr, int bondExpr, boolean quiet)
   {  super(oMolType, rmsdTag, atomDevTag, atomInfoTags, doOptimize, doMirror, quiet);

      int mcsstype = OEMCSType.Exhaustive;
      for( OEMol rMol : refmols)
      {  OEMCSSearch mcss = new OEMCSSearch(rMol, atomExpr, bondExpr, mcsstype);
         mcss.SetMCSFunc(new OEMCSMaxBondsCompleteCycles());
         mcsss.add(mcss);
      }
   }

   @Override
   public List<OEMolBase> getReference()
   {  ArrayList<OEMolBase> pats = new ArrayList<OEMolBase>(mcsss.size());
      for( OEMCSSearch ss : mcsss)
         pats.add(ss.GetPattern());

      return pats;
   }

   @Override
   public List<OEMatchBaseIter> getMatchIterator(OEMolBase fitMol)
   {  ArrayList<OEMatchBaseIter> its = new ArrayList<OEMatchBaseIter>(mcsss.size());
      for( OEMCSSearch ss : mcsss)
         its.add(ss.Match(fitMol,false));

      return its;
   }

   @Override
   public void close()
   {  for( OEMCSSearch ss : mcsss)
         ss.delete();
   }
}
