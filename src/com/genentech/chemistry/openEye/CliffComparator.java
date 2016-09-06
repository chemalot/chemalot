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

import java.util.regex.Pattern;

import openeye.oechem.OEGraphMol;
import openeye.oechem.OEMolBase;
import openeye.oechem.oechem;

import com.aestel.chemistry.openEye.SimComparator;
import com.genentech.chemistry.openEye.AAPathComparatorFact.AAPathCompareType;

/**
 * This comparator comapres two molecule objects and computes the activity cliff
 * of the two. The cliff is defined as  abs(activity1 - activity2)/(1-similarity.
 *
 * Activity is defined by an sdf tag value, the tagname being provided to the constructor.
 * The similarity is computed with the SimCoparator passed to the constructor.
 *
 * @author albertgo
 *
 */
public class CliffComparator implements SimComparator<OEMolBase>
{  private static final Pattern REMOVEModPat = Pattern.compile("^[<> ~=]+");
   private OEMolBase mol;
   private SimComparator<OEMolBase> simCompartor;
   private final boolean takeLog;
   private final String propertyTag;
   private final double minCliffSim;
   private final boolean ignoreModifier;
   final Double propertyValue;


   /**
    *
    * @param mol the molecule object is owned adn deleted by this comparator.
    */
   CliffComparator(SimComparator<OEMolBase> simComparator, OEMolBase mol, String propertyTag, boolean ignoreModifier, boolean takeLog, double minCliffSim)
   {  this.mol = mol;
      this.simCompartor = simComparator;
      this.takeLog = takeLog;
      this.ignoreModifier = ignoreModifier;
      this.propertyTag = propertyTag;
      this.minCliffSim = minCliffSim;
      propertyValue = getValue(mol);
   }

   @Override
   public double similarity(OEMolBase other)
   {  if( propertyValue == null ) return 0D;

      Double otherVal = getValue(other);
      if( otherVal == null ) return 0D;

      double sim = simCompartor.similarity(other);
      if( sim < minCliffSim ) return 0D;
      if( sim > 0.9999 )
         if ( Double.compare(propertyValue,otherVal) == 0 )
            return 0D;       // flag for identical compound
         else
            return 999D;       // max Cliff

      sim = computeCliff(otherVal, sim);
      if( sim >= 999D ) return 999D;
      return sim;
   }


   @Override
   public double similarity(SimComparator<OEMolBase> otherC)
   {  if( propertyValue == null ) return 0;

      CliffComparator other = (CliffComparator)otherC;
      Double otherVal = other.propertyValue;
      if( otherVal == null ) return 0;

      double sim = simCompartor.similarity(other.simCompartor);
      if( sim < minCliffSim ) return 0D;
      if( sim > 0.9999 )
         if ( Double.compare(propertyValue,otherVal) == 0 )
            return 0;       // flag for identical compound
         else
            return 999D;       // max Cliff


      sim = computeCliff(otherVal, sim);
      if( sim >= 999D ) return 999D;
      return sim;
   }

   /**
    * Default implementation to be improved.
    */
   @Override
   public double similarity(SimComparator<OEMolBase> other, double minSim)
   { return similarity(other);
   }


   double computeCliff(Double otherVal, double sim)
   {  return Math.abs(propertyValue - otherVal) / (1D - sim);
   }

   /**
    * This will also close the parent comparator.
    */
   @Override
   public void close()
   {  simCompartor.close();
      // there might be no need to delete mol as it will might have been deleted
      // by the simComparator
      if( mol != null && mol.IsValid()) mol.delete();
   }


   private Double getValue(OEMolBase mol2)
   {  String pval = oechem.OEGetSDData(mol2, propertyTag);
      if( ignoreModifier ) pval = REMOVEModPat.matcher(pval).replaceFirst("");

      if( pval == null || pval.length() == 0) return null;

      if( takeLog )
         return Math.log10(Double.parseDouble(pval));
      else
         return Double.parseDouble(pval);
   }

   public static void main(String argv[])
   {  String smiPat;
      String smiTar;

      // aa paper 1
      smiPat = "c1ccn[nH]1";
      smiTar = "Cc1ccno1";

      OEGraphMol mol1 = new OEGraphMol();
      oechem.OEParseSmiles(mol1, smiPat);
      oechem.OESetSDData(mol1, "prop", "1000");
      System.err.println(mol1.IsValid() + " " + mol1.NumAtoms());

      OEGraphMol mol2 = new OEGraphMol();
      oechem.OEParseSmiles(mol2, smiTar);
      oechem.OESetSDData(mol2, "prop", "10");
      System.err.println(mol2.IsValid() + " " + mol2.NumAtoms());

      AAPathComparatorFact aaPathFact = new AAPathComparatorFact(AAPathCompareType.DEFAULT,
                                                                 AAPathComparatorFact.DEFAULTVersion);
      OEGraphMol mol11 = new OEGraphMol(mol1);
      SimComparator<OEMolBase> simer = aaPathFact.createComparator(mol11);
      double sim   = simer.similarity(mol2);

      CliffComparatorFact cliffCFact = new CliffComparatorFact(aaPathFact, "prop", false, true, 0.3D, 0, 0);
      mol11 = new OEGraphMol(mol1);
      CliffComparator cliffer = cliffCFact.createComparator(mol11);
      double cliff = cliffer.similarity(mol2);
      System.out.printf("sim=%g delta=%d cliff=%f\n", sim, 2, cliff);
      cliffer.close();
      cliffCFact.close();

      cliffCFact = new CliffComparatorFact(aaPathFact, "prop", false, true, 0.3D, 100, 0.5);
      mol11 = new OEGraphMol(mol1);
      cliffer = cliffCFact.createComparator(mol11);
      cliff = cliffer.similarity(mol2);
      System.out.printf("foldAt1Sim=%f simAt2Delta=%f sim=%g delta=%d cliff=%f\n", 100D, .5, sim, 2, cliff);
      cliffer.close();
      cliffCFact.close();

      cliffCFact = new CliffComparatorFact(aaPathFact, "prop", false, true, 0.3D, 100, 0.4);
      mol11 = new OEGraphMol(mol1);
      cliffer = cliffCFact.createComparator(mol11);
      cliff = cliffer.similarity(mol2);
      System.out.printf("foldAt1Sim=%f simAt2Delta=%f sim=%g delta=%d cliff=%f\n", 100D, .4, sim, 2, cliff);
      cliffer.close();
      cliffCFact.close();


      mol11.delete();
      mol2.delete();
   }
}
