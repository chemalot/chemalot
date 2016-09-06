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

import openeye.oechem.OEAtomBase;
import openeye.oechem.OEAtomBaseIter;
import openeye.oechem.OEAtomBaseVector;
import openeye.oechem.OEAtomStereo;
import openeye.oechem.OEBondBase;
import openeye.oechem.OEBondBaseIter;
import openeye.oechem.OEBondStereo;
import openeye.oechem.OEExprOpts;
import openeye.oechem.OEGraphMol;
import openeye.oechem.OEMolBase;
import openeye.oechem.OEQMol;
import openeye.oechem.OESMILESFlag;
import openeye.oechem.OESubSearch;
import openeye.oechem.OEUnaryAtomPred;
import openeye.oechem.OEUnaryBondPred;
import openeye.oechem.oechem;

import com.genentech.oechem.tools.Atom;
import com.genentech.oechem.tools.AtomHasStereoSpecifiedFunctor;

/**
 * This class provides a function to randomly reorder atoms in an OEGrapMol.
 *
 * This is to be used to validate the smiles canonicalization routine in oechem.
 * @author albertgo
 *
 */
public class MyCannonicalize {
   private static final int NCHECK = 20;

   private final OEGraphMol tmpMol;
   private final OEQMol tmpQMol;
   private final OEUnaryAtomPred atHasStereoSpecFunct;
   private final OEUnaryBondPred bdHasStereoSpecFunc;
   private final AtomStereoNonChiralFunctor atNonStereoChiralFunc;
   private final OEAtomBaseVector tmpAtVec;

   public MyCannonicalize() {
      tmpMol = new OEGraphMol();
      atHasStereoSpecFunct = new AtomHasStereoSpecifiedFunctor();
      atNonStereoChiralFunc= new AtomStereoNonChiralFunctor();
      bdHasStereoSpecFunc  = new BondHasStereoSpecifiedFunctor();
      tmpQMol = new OEQMol();
      tmpAtVec = new OEAtomBaseVector();
   }


   /**
    * Create canonical Smiles adding additional checks to the openeye functionality
    * in order to fix shortcomings:
    *
    * - If mol has stereo atoms invert stereo centers and try atom atom match, if
    *   that works return lexically smaller smiles to circumvent OEBug on eg.
    *   1,4 Dimethyl cyclohexan.
    * - Invert all stereo centers non no-chiral atoms and try atom atom match if
    *   that works return lexically smaller smiles.
    * - invert all double bonds if atom atom match yields same compound return
    *   lexically smaller smiles.
    * - generate NCHECK molecule objects with random atom and bond orderings,
    *   compare for smallest smiles.
    *
    * @param isoSmi if true an isomeric smiles will be created.
    */
   public String myCanonicalSmi(OEMolBase mol, boolean isoSmi) {
      oechem.OESuppressHydrogens(mol,false,false,true);
      oechem.OEPerceiveChiral(mol);

      int smiFlavor = OESMILESFlag.Canonical;
      if( isoSmi ) smiFlavor |= OESMILESFlag.AtomStereo | OESMILESFlag.BondStereo
                            | OESMILESFlag.Isotopes;
      String minSmi = oechem.OECreateSmiString(mol, smiFlavor);
      tmpMol.Clear();

      boolean weTrustOEChemStereoSmilesGeneration = true;
      if( ! weTrustOEChemStereoSmilesGeneration && isoSmi) {
         oechem.OEAddMols(tmpQMol, mol);

         tmpQMol.BuildExpressions(OEExprOpts.ExactAtoms, OEExprOpts.ExactBonds);
         tmpQMol.BuildExpressions(OEExprOpts.SymmetryClass, OEExprOpts.ExactBonds);
         tmpQMol.BuildExpressions(OEExprOpts.SymmetryClass, OEExprOpts.SymmetryClass);
         tmpQMol.BuildExpressions(OEExprOpts.ExactAtoms, OEExprOpts.ExactBonds);
         OESubSearch subSearch = new OESubSearch(tmpQMol);

         // invert all atoms with stereochemistry and check for lowest smiles
         // catches C[C@H]1CCC[C@H]1C
         if( invertAtoms(atHasStereoSpecFunct, tmpMol) > 1) {
            if(subSearch.SingleMatch(tmpMol))
            {  // this is a smiles which did not get canonized correctly
               String smi2 = oechem.OECreateSmiString(tmpMol, smiFlavor);
               if(smi2.compareTo(minSmi) < 0 ) {
                  System.err.printf("A Unstable Canonicalization %s>%s\n", minSmi,smi2);
                  minSmi = smi2;
               }
            }
         }

         // invert all stereo centers which have stereo but are not chiral
         // catches C[C@H]1CC[C@H](CC1)CC[C@@H](C)Cl
         tmpMol.Clear();
         oechem.OEAddMols(tmpMol, mol);
         if( invertAtoms(atNonStereoChiralFunc, tmpMol) > 1) {
            if(subSearch.SingleMatch(tmpMol))
            {  // this is a smiles which did not get canonized correctly
               String smi2 = oechem.OECreateSmiString(tmpMol, smiFlavor);
               if(smi2.compareTo(minSmi) < 0 ) {
                  System.err.printf("B Unstable Canonicalization %s>%s\n", minSmi,smi2);
                  minSmi = smi2;
               }
            }
         }

         // invert all double bonds
         // catches C/C(=N/OC(=O)CCC(=O)O/N=C(\C)/c1cccs1)/c2cccs2
         tmpMol.Clear();
         oechem.OEAddMols(tmpMol, mol);
         if( invertBonds(bdHasStereoSpecFunc, tmpMol) > 0) {
            if(subSearch.SingleMatch(tmpMol))
            {  // this is a smiles which did not get canonized correctly
               String smi2 = oechem.OECreateSmiString(tmpMol, smiFlavor);
               if(smi2.compareTo(minSmi) < 0 ) {
                  System.err.printf("C Unstable Canonicalization %s>%s\n", minSmi,smi2);
                  minSmi = smi2;
               }
            }
         }
         subSearch.delete();

         tmpMol.Clear();
      }

      oechem.OEAddMols(tmpMol, mol);
      // now generate NCHECK mol objects with random ordering and check for
      // smallest smiles
      for(int i=1; i<=NCHECK; i++) {
         oechem.OEScrambleMolecule(tmpMol);

         String smi2 = oechem.OECreateSmiString(tmpMol, smiFlavor);
         if( smi2.compareTo(minSmi) < 0 ) {
            System.err.printf("D Unstable Canonicalization %s>%s\n", minSmi,smi2);
            minSmi = smi2;
         }
      }

      return minSmi;
   }

   /** invalidates this instance of {@link MyCannonicalize}
    */
   public void delete() {
      atHasStereoSpecFunct.delete();
      atNonStereoChiralFunc.delete();
      bdHasStereoSpecFunc.delete();
      tmpMol.delete();
      tmpQMol.delete();
      tmpAtVec.delete();
   }


   /** Invert double bond stereochemistry on bonds in bondFunctor */
   private int invertBonds(OEUnaryBondPred bondFunctor, OEMolBase mol) {
      int count =0;
      OEBondBaseIter bIt = mol.GetBonds(bondFunctor);
      while( bIt.hasNext() ) {
         OEBondBase bd = bIt.next();
         if(! bd.HasStereoSpecified() ) continue;
         if( bd.IsInRing() ) continue;

         OEAtomBase at1 = bd.GetBgn();
         OEAtomBase at2 = bd.GetEnd();
         tmpAtVec.clear();
         tmpAtVec.add(Atom.findExoNeighbor(at1, at2));
         tmpAtVec.add(Atom.findExoNeighbor(at2, at1));
         int stereo = bd.GetStereo(tmpAtVec, OEBondStereo.CisTrans);
         stereo = stereo == OEBondStereo.Cis ? OEBondStereo.Trans : OEBondStereo.Cis;
         bd.SetStereo(tmpAtVec, OEBondStereo.CisTrans, stereo);

         count++;
      }
      bIt.delete();

      return count;
   }


   /** Invert atoms identified by atomFunctor */
   private int invertAtoms(OEUnaryAtomPred atomFunctor, OEMolBase mol) {
      int count =0;
      OEAtomBaseIter aIt = mol.GetAtoms(atomFunctor);
      while( aIt.hasNext() ) {
         OEAtomBase at = aIt.next();
         invertStereo(at);
         count++;
      }
      aIt.delete();

      return count;
   }


   private void invertStereo(OEAtomBase at) {
      tmpAtVec.clear();

      OEAtomBaseIter neighborsIt=at.GetAtoms();
      while( neighborsIt.hasNext() ) {
         OEAtomBase atNeighbor = neighborsIt.next();
         tmpAtVec.add(atNeighbor);
      }
      neighborsIt.delete();

      int stereo = at.GetStereo(tmpAtVec, OEAtomStereo.Tetrahedral);
      if(stereo == OEAtomStereo.Left )
         stereo = OEAtomStereo.Right;
      else
         stereo = OEAtomStereo.Left;

      at.SetStereo(tmpAtVec, OEAtomStereo.Tetrahedral, stereo);
   }

   public static final void main(String...args)
   {  OEGraphMol mol = new OEGraphMol();
      String smi = "C[C@@H]2[C@@H]([C@H]([C@H]2Cl)C)Cl";
      smi = "C[C@H]1CCC[C@H]1C";

      oechem.OEParseSmiles(mol, smi);
      String canSmi = smi;

      MyCannonicalize canonizer = new MyCannonicalize();
      for(int i=0; i<1000; i++) {
         oechem.OEScrambleMolecule(mol);

         String smiOut = canonizer.myCanonicalSmi(mol,true);
         if(canSmi.equals(smiOut)) continue;

         System.err.printf("%4d %s\n",i,smiOut);
         canSmi = smiOut;
      }
      canonizer.delete();
   }
}
