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

import java.util.List;

import com.aestel.utility.DataFormat;
import com.genentech.chemistry.tool.align.SDFAlign.OUTType;

import openeye.oechem.*;

class FSSAlign implements AlignInterface
{  final List<OEMol> refMols;
   protected final boolean doMirror;
   protected final String rmsdTag;
   protected final boolean doOptimize;

   public FSSAlign(List<OEMol> refmol, OUTType oMolType, String rmsdTag, boolean doOptimize, boolean doMirror)
   {  this.doMirror = doMirror;
      this.rmsdTag = rmsdTag;
      this.doOptimize = doOptimize;
      this.refMols = refmol;
   }


   @Override
   public void align(OEMolBase fitmol)
   {  double rmat[] = new double[9];
      double trans[] = new double[3];
      double rmat2[] = new double[9];
      double trans2[] = new double[3];
      double minRMSD = 9999;
      OEGraphMol mirMol = null;
      boolean isMirror = false;

      oechem.OEPerceiveChiral(fitmol);
      if( doMirror && ! AbstractAlign.isChiral(fitmol) )
      {  mirMol = new OEGraphMol(fitmol);
         AbstractAlign.createMirror(mirMol);
      }

      for( OEMolBase refMol: refMols)
      {  double newRMSD = oechem.OERMSD(refMol, fitmol, true, true, doOptimize, rmat2, trans2);
         if( newRMSD < 0 ) // no match
            continue;

         if( newRMSD < minRMSD )
         {  minRMSD = newRMSD;
            rmat = rmat2;
            trans = trans2;
            isMirror = false;
         }

         if( mirMol != null )
         {  newRMSD = oechem.OERMSD(refMol, mirMol, true, true, doOptimize, rmat2, trans2);

            if( newRMSD < minRMSD )
            {  minRMSD = newRMSD;
               rmat = rmat2;
               trans = trans2;
               isMirror = true;
            }
         }
      }

      if( minRMSD < 9999 ) // no match
         return;

      if( isMirror )
      {  fitmol.Clear();
         oechem.OEAddMols(fitmol, mirMol);
         fitmol.SetTitle(mirMol.GetTitle());
      }
      if( mirMol != null) mirMol.delete();

      if( doOptimize )
      {  oechem.OERotate(fitmol, rmat);
         oechem.OETranslate(fitmol, trans);
      }
      if( rmsdTag != null )
          oechem.OESetSDData(fitmol, rmsdTag, DataFormat.formatNumber(minRMSD, "si3"));
   }


   @Override
   public void close()
   {  for( OEMolBase refMol: refMols)
         refMol.delete();
   }

}
