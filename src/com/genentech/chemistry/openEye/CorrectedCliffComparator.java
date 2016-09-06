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

import com.aestel.chemistry.openEye.SimComparator;

/**
/**
 * This comparator comapres two molecule objects and computes the activity cliff
 * of the two.
 * the Cliff calculation is corrected to de-emphasize low similarities and
 * emphasize medium similarities.
 *
 * Compute cliffs as abs(delta(prop)) * e^(a*sim^b) where:
 * a=ln(foldAt1Sim)
 * b=log(ln(2)/a), simAt2Delta)
 * simAt2Delta = similarity at which delta(prop) is doubled.
 * foldAt1Sim  = Maximum Cliff size is foldAt1Sim * delta(props)
 *
 * cf. CliffCorrection.xlsx
 *
 * @author albertgo
 *
 */
public class CorrectedCliffComparator extends CliffComparator
{  private final double a;
   private final double b;

   CorrectedCliffComparator(SimComparator<OEMolBase> simComparator, OEMolBase mol,
            String propertyTag, boolean ignoreModifier, boolean takeLog, double minCliffSim,
            double foldAt1Sim, double simAt2Delta)
   {  super(simComparator, mol, propertyTag, ignoreModifier, takeLog, minCliffSim);
      this.a = Math.log(foldAt1Sim);
      this.b = Math.log(Math.log(2)/a)/Math.log(simAt2Delta);
   }

   @Override
   double computeCliff(Double otherVal, double sim)
   {  return Math.abs(propertyValue - otherVal) * Math.exp(a*Math.pow(sim, b));
   }
}
