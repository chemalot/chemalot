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

import openeye.oechem.*;

public class LipophilicTyper implements AtomTyperInterface
{  private final int idx;

   public LipophilicTyper(int idx)
   {  this.idx =idx;
   }


   @Override
   public boolean isType(OEAtomBase at)
   {  switch(at.GetAtomicNum())
      {  case OEElemNo.C:
         case OEElemNo.Cl:
            // only C or CL can by Lipophilic and only if they have less than 4 non H neighbors
            if( at.GetExplicitDegree() >= 4 ) return false;


            OEAtomBaseIter atIt = at.GetAtoms();
            boolean isLipo = true;
            while(atIt.hasNext())
            {  OEAtomBase at2 = atIt.next();
               switch(at2.GetAtomicNum())
               {  case OEElemNo.C:
                     OEAtomBaseIter atIt2 = at2.GetAtoms();
                     // May not have polar neighbors (only C,F,,LC,Br are allowed
                     while(atIt2.hasNext())
                     {  OEAtomBase at3 = atIt2.next();
                        switch( at3.GetAtomicNum() )
                        {  case OEElemNo.C:
                           case OEElemNo.F:
                           case OEElemNo.Cl:
                           case OEElemNo.Br:
                              break;
                           default:
                              // if conjugated via double or aromatic bond to heteroatom
                              if( at2.GetBond(at3).GetOrder() % 2 == 0
                                 || at.GetBond(at2).GetOrder() % 2 == 0 )
                                 isLipo = false;
                           break;
                        }
                        if( ! isLipo ) break;
                     }
                     atIt2.delete();
                     break;

                  case OEElemNo.Cl:
                  case OEElemNo.F:
                  case OEElemNo.Br:
                     break;
                  default:
                     isLipo = false;
                     break;
               }
               if( ! isLipo ) break;
            }
            atIt.delete();
            return isLipo;

         default:
            return false;
      }
   }


   @Override
   public int getTypeIdx()
   {  return idx;
   }


   @Override
   public String getTypeName()
   {  return "L";
   }


   @Override
   public void close()
   {  // nothing to do
   }
}
