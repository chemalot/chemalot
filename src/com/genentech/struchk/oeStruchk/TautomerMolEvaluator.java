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

import openeye.oechem.*;

import com.genentech.oechem.tools.*;


public class TautomerMolEvaluator {
   private static final int NONCHIRALNotSpecified  = 0;   // must be 0 because this is default
   private static final int CHIRALSpecified        = 1;
   private static final int CHIRALNotSpecified     = 2;
   private static final int NONCHIRALSpecified     = 3;

   private static final int DUMMYTag = oechem.OEGetTag("dummy");

   private final OEMolBase bestMol = new OEGraphMol();
   private final OEMolBase tmpMol  = new OEGraphMol();
   private final OEUnaryBondPred dBondFunct = new BondOrderFunctor(2);
   private final OEUnaryBondPred chiralBondFunct = new BondIsChiralFunctor();
   private final OEUnaryBondPred chiralDBondFunct = new OEAndBond(dBondFunct, chiralBondFunct);
   private final OEUnaryAtomPred atomHasStereoFunct = new AtomHasStereoSpecifiedFunctor();
   private final OEUnaryAtomPred atomIsChiralFunct = new AtomIsChiralFunctor();
   private final OEUnaryAtomPred stereoAtomFunct = new OEOrAtom(atomIsChiralFunct, atomHasStereoFunct);

   private int bestAromaticCount = -1;
   private int bestChargeCount = 0;
   private String bestSmi;
   private int molCount;
   private volatile boolean isDeleted = false;

   public TautomerMolEvaluator() {
   }

   /**
    * Be sure to call reset before using TautomerMolFunctor on a new molecule.
    */
   public boolean evaluate(OEMolBase in)
   {  //System.err.println("evaluate: "+ oechem.OECreateIsoSmiString(in));

      tmpMol.Clear();
      oechem.OEAddMols(tmpMol, in);
      oechem.OEPerceiveChiral(tmpMol);
      //oechem.OEAssignMDLHydrogens(tmpMol);

      if( checkForStereoChanges(tmpMol))
      {  // System.err.printf("Stereo Changed: %s\n", oechem.OECreateSmiString(tmpMol),OESMILESFlag.ISOMERIC);
         return false;
      }

      int aromaticCount = 0;
      int chargesCount = 0;
      oechem.OEAssignAromaticFlags(tmpMol);
      OEAtomBaseIter aIt = tmpMol.GetAtoms();
      while( aIt.hasNext() ) {
         OEAtomBase at = aIt.next();
         chargesCount += Math.abs(at.GetFormalCharge());

         if( at.IsAromatic() )
         {  int atomScore = 0;
            OEBondBaseIter bIt = at.GetBonds();
            while ( bIt.hasNext() ) {
               if( bIt.next().IsAromatic() ) atomScore++;
            }
            bIt.delete();
            aromaticCount += atomScore-1;
         }
      }
      aIt.delete();
      if( chargesCount < bestChargeCount
            || (chargesCount == bestChargeCount && aromaticCount >= bestAromaticCount)
            || molCount == 0 ) {
         String newSmi = oechem.OECreateIsoSmiString(tmpMol);
         // System.err.printf("newS=%s\tbestS=%s\n",newSmi, bestSmi);
         if( molCount == 0
             || chargesCount < bestChargeCount
             || aromaticCount > bestAromaticCount
             || newSmi.compareTo(bestSmi) > 0 ) {
            bestSmi = newSmi;
            bestMol.Clear();
            oechem.OEAddMols(bestMol, in);
            // System.err.printf("Keeping mol: %s\n", OETools.molToCanSmi(bestMol, true));
            bestAromaticCount  = aromaticCount;
            bestChargeCount = chargesCount;
            molCount++;
         }
      }

      // System.err.println("call: "+ (tautCount < maxStates) + " " + oechem.OECreateIsoSmiString(in));
      // printAtoms(in);
      return true;
   }

   private boolean checkForStereoChanges(OEMolBase in) {
      boolean stereoChanged = false;

      oechem.OEAssignAromaticFlags(in);
      oechem.OEPerceiveChiral(in);
      
      // Disregard tautomers changing double bond stereochemistry:
      // C1(C2=NNC=C2)=CC=CN1>>C3(/C=CC=N3)=C4C=CNN/4
      
      OEBondBaseIter bdIt = in.GetBonds();
      while( bdIt.hasNext() ) {
          OEBondBase bd = bdIt.next();

          int bdStereoType  = getStereoType(bd);
          int orgStereoType = bd.GetIntData(OEStruchk.HADStereoTag);
          if( bdStereoType != orgStereoType ) {
              //System.err.printf("new=%s orig=%s bond=%s ism=%s\n",
              //        stereoFragToString(bdStereoType), stereoFragToString(orgStereoType),
              //        Bond.toString(bd), oechem.OEMolToSmiles(in));
              stereoChanged = true;
              break;
          }
      }
      bdIt.delete();
      if( stereoChanged ) {
         return true;
      }
      
      // capture cases like this:
      // CC1=CC(C(C(C2=C(O)OC(C)=CC2=O)C)=C(O)O1)=O>>CC3=CC(O)=C(C(C4=C(O)C=C(C)OC4=O)C)C(O3)=O
      // the central atom is chiral on the left but not on the right

      // flag atoms which had stereo on input
      OEAtomBaseIter aIt = in.GetAtoms();
      while( aIt.hasNext() ) {
         OEAtomBase at = aIt.next();

         int atStereoType  = getStereoType(at);
         int orgStereoType = at.GetIntData(OEStruchk.HADStereoTag);
         if( atStereoType != orgStereoType ) {
            //System.err.printf("new=%s orig=%s atom=%s ism=%s\n",
            //        stereoFragToString(atStereoType), stereoFragToString(orgStereoType),
            //        Atom.toString(at), oechem.OEMolToSmiles(in));
            stereoChanged = true;
            break;
         }
      }
      aIt.delete();

      return stereoChanged;
   }

  /**
    * Setup flags for atoms and bonds with stereo so that we can recognize changes.
    */
   public void setupNewMol(OEMolBase in) {
      oechem.OEPerceiveChiral(in);
      bestSmi = null;
      bestMol.Clear();
      oechem.OEAddMols(bestMol, in);
      bestAromaticCount = -1;
      bestChargeCount = 0;
      molCount = 0;

      // flag bonds which had stereo on input
      OEBondBaseIter bdIt = in.GetBonds(chiralDBondFunct);
      while( bdIt.hasNext() ) {
         OEBondBase bd = bdIt.next();

         int bdStereoType = getStereoType(bd);
         if( bdStereoType != NONCHIRALNotSpecified ) {
            bd.SetIntData(OEStruchk.HADStereoTag, bdStereoType);
            // System.err.printf("has bond stereo %s %s \n", Bond.toString(bd), stereoFragToString(bdStereoType));
         }
      }
      bdIt.delete();

      // flag atoms which had stereo on input
      OEAtomBaseIter aIt = in.GetAtoms(stereoAtomFunct);
      while( aIt.hasNext() ) {
         OEAtomBase at = aIt.next();

         int atStereoType = getStereoType(at);
         if( atStereoType != NONCHIRALNotSpecified ) {
            at.SetIntData(OEStruchk.HADStereoTag, atStereoType);
            // System.err.printf("has atom stereo %s %s \n", Atom.toString(at), stereoFragToString(atStereoType));
         }
      }
      aIt.delete();
   }

   private static int getStereoType(OEBondBase bd) {
       if (bd.GetBoolData(OEStruchk.STEREOClearTag)) {
           assert (! bd.HasStereoSpecified())
           : "Bond stereo should have been removed!";
           return 0;
       }

       int bdStereoType;
       if (bd.IsChiral())
           if (bd.HasStereoSpecified(OEBondStereo.CisTrans))
               bdStereoType = CHIRALSpecified;
           else
               bdStereoType = CHIRALNotSpecified;
       else if (bd.HasStereoSpecified(OEBondStereo.CisTrans))
           bdStereoType = NONCHIRALSpecified;
       else
           bdStereoType = NONCHIRALNotSpecified;
       return bdStereoType;
   }
   
   private static int getStereoType(OEAtomBase at) {
      if( at.GetBoolData(OEStruchk.STEREOClearTag) ) {
         assert (! at.HasStereoSpecified()) || at.GetBoolData(OEStruchk.ATROPIsomericCenter)
               : "Atom stereo should have been removed!";
         return 0;
      }

      int atStereoType;
      if( at.IsChiral() )
         if( at.HasStereoSpecified(OEAtomStereo.Tetrahedral))
            atStereoType = CHIRALSpecified;
         else
            atStereoType = CHIRALNotSpecified;
      else
         if( at.HasStereoSpecified(OEAtomStereo.Tetrahedral))
            atStereoType = NONCHIRALSpecified;
         else
            atStereoType = NONCHIRALNotSpecified;
      return atStereoType;
   }

   /**
    * @return a reference to an internal molecule do not delete or modify
    */
   public OEMolBase GetBest()
   {  return bestMol;
   }

   @Override
   public void finalize()
   {  delete();
   }

   public void delete()
   {  if(isDeleted) return;
      isDeleted  = true;

      bestMol.delete();
      tmpMol.delete();
      chiralDBondFunct.delete();
      dBondFunct.delete();
      chiralBondFunct.delete();
      stereoAtomFunct.delete();
      atomHasStereoFunct.delete();
      atomIsChiralFunct.delete();
   }

   public static final String stereoFragToString(int stereoflag) {  
       switch( stereoflag ) {
           case NONCHIRALNotSpecified : return "NonChiralNoSpec";
           case CHIRALSpecified       : return "ChiralSpec";
           case CHIRALNotSpecified    : return "ChiralNotSpec";
           case NONCHIRALSpecified    : return "NonChiralNotSpec";
       }
       return "Unknown";
   }
   
   final void printAtoms(OEMolBase mol)
   {  OEAtomBaseIter aIt = mol.GetAtoms();
      System.err.print("Atom:\t");
      while( aIt.hasNext() ) {
         OEAtomBase at = aIt.next();
         System.err.printf("%s %d\t",Atom.getAtomName(at), at.GetIntData(OEStruchk.HADStereoTag));
      }
      aIt.delete();
      OEBondBaseIter bdIt = mol.GetBonds();
      System.err.print("Bond:\t");
      while( bdIt.hasNext() ) {
         OEBondBase bd = bdIt.next();
         System.err.printf("%s %b\t", bd.GetOrder(), bd.GetBoolData(OEStruchk.HADStereoTag));
      }
      bdIt.delete();
      System.err.println(mol.GetBoolData(DUMMYTag)+"\t"+oechem.OECreateIsoSmiString(mol));
   }
}
