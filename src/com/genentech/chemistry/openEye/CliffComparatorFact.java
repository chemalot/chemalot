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
import com.aestel.chemistry.openEye.SimComparatorFactory;

public class CliffComparatorFact implements SimComparatorFactory<OEMolBase, OEMolBase, SimComparator<OEMolBase>>
{  private final SimComparatorFactory<OEMolBase, OEMolBase, SimComparator<OEMolBase>> parentFact;
   private final boolean takeLog;
   private final String propertyTag;
   private final double minCliffSim;
   private final double foldAt1Sim;
   private final double simAt2Delta;
   private final boolean ignoreModifier;

   /**
    * Create a cliff comparator using the SimilarityComaprator provided.
    * The cliff value is calculated as delta(property)/(1-Sim).
    * If takelog is given formula is delta(log(property)/1-Sim)
    * @param parentFact Factory to create a new SimilarityComparators for the underlying
    *                   Similarity.
    * @param propertyTag Tag of property to use for numerator
    * @param takeLog if true the log of the property value is used
    * @param minCliffSim if Sim < minCliffSim return 0
    * @param foldAt1Sim  if not 0 use modified cliff formula cf. {@link CorrectedCliffComparator}.
    * @param simAt2Delta if not 0 use modified cliff formula cf. {@link CorrectedCliffComparator}.
    */
   public CliffComparatorFact(SimComparatorFactory<OEMolBase, OEMolBase, SimComparator<OEMolBase>> parentFact,
                              String propertyTag, boolean ignoreModifier, boolean takeLog, double minCliffSim,
                              double foldAt1Sim, double simAt2Delta)
   {  this.parentFact = parentFact;
      this.propertyTag = propertyTag;
      this.ignoreModifier = ignoreModifier;
      this.takeLog = takeLog;
      this.minCliffSim = minCliffSim;
      this.foldAt1Sim = foldAt1Sim;
      this.simAt2Delta = simAt2Delta;
   }

   /** returns new objects which should be deleted separately */
   @Override
   public OEMolBase createComparable(OEMolBase in)
   {  return parentFact.createComparable(in);
   }

   @Override
   public CliffComparator createComparator(OEMolBase mol)
   {  SimComparator<OEMolBase> parentComparator = parentFact.createComparator(mol);
      if( foldAt1Sim <= 0D || simAt2Delta <= 0D)
         return new CliffComparator(parentComparator, mol, propertyTag, ignoreModifier, takeLog, minCliffSim );

      return new CorrectedCliffComparator(parentComparator, mol, propertyTag, ignoreModifier, takeLog, minCliffSim, foldAt1Sim, simAt2Delta);
   }

   @Override
   public void close()
   {  parentFact.close();
   }
}
