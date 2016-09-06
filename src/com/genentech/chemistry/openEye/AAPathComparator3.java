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

import openeye.oechem.OEMolBase;

/**
 * Modified algorithm to compute AtomAtomPath similarity.
 *
 * The algorithm to compute the similarity of two atoms is changed as follows:
 *    pathsInCommon / max(nPath(AtomMolA), nPath(AtomMolB))
 *
 * This should prevent the situation where a small substituent results in a higher
 * similarity than a big substituent.
 *
 * @author albertgo
 *
 */
public class AAPathComparator3 extends AAPathComparator
{  AAPathComparator3(OEMolBase mol, int maxBonds, double minAtSim)
   {  super(mol, maxBonds, minAtSim);
   }

   @Override
   double computeAtomSim(int nPaths, int nPaths2, int common)
   {  return common / (double)Math.max(nPaths, nPaths2);
   }


   @Override
   double computeMoleculeSimilarity(int nAtoms1, int nAtoms2, double simSum)
   {  return simSum / Math.max(nAtoms1, nAtoms2);
   }

}
