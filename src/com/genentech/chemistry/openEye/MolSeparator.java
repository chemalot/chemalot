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
package com.genentech.chemistry.openEye;

import java.util.ArrayList;
import java.util.List;

import com.genentech.oechem.tools.OETools;

import openeye.oechem.OEGraphMol;
import openeye.oechem.OEMolBase;
import openeye.oechem.OEPartPredAtom;
import openeye.oechem.oechem;


/**
 * Separate components in a molfile into individual molecules.
 *
 * @author Man-Ling Lee / December 20, 2012
 * Copyright 2012-2013 Genentech
 */
public class MolSeparator
{  int nComponents;
   List<OEMolBase> components = null;

   public MolSeparator()
   {
      nComponents = 0;
      components = new ArrayList<OEMolBase>();
   }


   public void clear()
   {
      for( int i = 1; i < nComponents; i++ )
         components.get(i).delete();

      nComponents = 0;
      components.clear();
   }


   public int separate( OEMolBase inMol )  {
      int[] parts = new int[ inMol.GetMaxAtomIdx() ];
      nComponents = oechem.OEDetermineComponents( inMol, parts );

      OEPartPredAtom pred = new OEPartPredAtom( parts );

      // get components from input molecule that can consists of one or more
      // components
      for( int i = 1; i <= nComponents; i++ )
      {
        pred.SelectPart(i);
        OEGraphMol partMol = new OEGraphMol();
        oechem.OESubsetMol( partMol, inMol, pred );
        components.add( partMol );
      }
      pred.delete();

      return nComponents;
   }


   public OEMolBase getMol( int index )
   {
      if( index > components.size() )
         throw new IndexOutOfBoundsException(
                     "Input molfile contains only " + nComponents + " components" );
      return components.get( index );
   }


   public String getCanISmi( int index )
   {
      if( index > components.size() )
         throw new IndexOutOfBoundsException(
                     "Input molfile contains only " + nComponents + " components" );
      return OETools.molToCanSmi( components.get( index ), true );
   }
}
