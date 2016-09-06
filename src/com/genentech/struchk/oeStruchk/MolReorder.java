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

import java.util.Random;

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
public class MolReorder {
   private static final Random rand = new Random();
   private static final int NCHECK = 20;
   private static int outFlavor = OESMILESFlag.AtomStereo | OESMILESFlag.BondStereo
                           | OESMILESFlag.Isotopes | OESMILESFlag.Hydrogens
                           | OESMILESFlag.Canonical;

   private MolReorder() {
   }

   public static OEGraphMol reorder(OEGraphMol inMol) {
      OEAtomBase[] ats = new OEAtomBase[inMol.NumAtoms()];

      int atIdx = 0;
      OEAtomBaseIter aIt = inMol.GetAtoms();
      while( aIt.hasNext() )
         ats[atIdx++] = aIt.next();
      aIt.delete();

      // generate random atom sequence
      for(int from=0; from<ats.length; from++) {
         int to   = rand.nextInt(ats.length);
         OEAtomBase at = ats[from];
         ats[from] = ats[to];
         ats[to]   = at;
      }

      // create new atoms and keep copy in the order of the old molecule
      OEGraphMol tmpMol = new OEGraphMol();
      OEAtomBase[] newAtomsAtOldPos = new OEAtomBase[inMol.GetMaxAtomIdx()];
      for(int i=0; i<ats.length; i++) {
         if( ats[i] == null) continue;

         newAtomsAtOldPos[ats[i].GetIdx()] = tmpMol.NewAtom(ats[i]);
      }

      // generate random bond sequence
      OEBondBase[] bds = new OEBondBase[inMol.NumBonds()];

      int bdIdx = 0;
      OEBondBaseIter bIt = inMol.GetBonds();
      while( bIt.hasNext() )
         bds[bdIdx++] = bIt.next();
      bIt.delete();

      // generate bond random sequence
      for(int from=0; from<bds.length; from++) {
         int to   = rand.nextInt(bds.length);
         OEBondBase bd = bds[from];
         bds[from] = bds[to];
         bds[to]   = bd;
      }

      // create new bonds and keep copy in order of old molecule
      OEBondBase[] newBondsInOldPos = new OEBondBase[inMol.GetMaxBondIdx()];
      for(int i=0; i<bds.length; i++ ) {
         OEBondBase bd = bds[i];
         if(bd == null) continue;

         OEAtomBase at1 = newAtomsAtOldPos[bd.GetBgn().GetIdx()];
         OEAtomBase at2 = newAtomsAtOldPos[bd.GetEnd().GetIdx()];
         if(rand.nextBoolean()) {   // random swap
            OEAtomBase tmpAt = at1;
            at1 = at2;
            at2 = tmpAt;
         }

         newBondsInOldPos[bd.GetIdx()] =
               tmpMol.NewBond(at1, at2, bd.GetOrder() );
      }

      // copy tetrahedral chirality from old to new
      OEAtomBaseVector oldVec = new OEAtomBaseVector();
      OEAtomBaseVector newVec = new OEAtomBaseVector();
      aIt = inMol.GetAtoms();
      while( aIt.hasNext() ) {
         OEAtomBase oldAt = aIt.next();
         if( ! oldAt.HasStereoSpecified(OEAtomStereo.Tetrahedral) ) continue;

         OEAtomBase newAt = newAtomsAtOldPos[oldAt.GetIdx()];

         for (OEAtomBaseIter neighborsIt=oldAt.GetAtoms(); neighborsIt.hasNext();) {
            OEAtomBase oldAtNeighbor = neighborsIt.next();
            OEAtomBase newAtNeighbor = newAtomsAtOldPos[oldAtNeighbor.GetIdx()];
            oldVec.add(oldAtNeighbor);
            newVec.add(newAtNeighbor);
         }
         int stereo = oldAt.GetStereo(oldVec, OEAtomStereo.Tetrahedral);
         newAt.SetStereo(newVec, OEAtomStereo.Tetrahedral, stereo);
         newVec.clear();
         oldVec.clear();
      }
      aIt.delete();

      // copy double bond stereochemistry from old to new
      bIt = inMol.GetBonds();
      while( bIt.hasNext() ) {
         OEBondBase oldBd = bIt.next();

         if(! oldBd.HasStereoSpecified(OEBondStereo.CisTrans)) continue;
         OEBondBase newBd = newBondsInOldPos[oldBd.GetIdx()];

         // SetStereo needs two atoms connected to the beginning and end point of
         // a double bond in a vector.
         OEAtomBase at1 = oldBd.GetBgn();
         OEAtomBase at2 = oldBd.GetEnd();

         // find atom on beginning atom
         OEAtomBase oldNat = Atom.findExoNeighbor(at1, at2);    // find exocyclic bond from at1, not being at2
         if( oldNat == null) continue;   // should not happen bond has no other bonds

         oldVec.add(oldNat);
         newVec.add(newAtomsAtOldPos[oldNat.GetIdx()]);

         // find atom on end atom
         oldNat = Atom.findExoNeighbor(at2, at1);    // find exocyclic bond from at1, not being at2
         if( oldNat == null) {
            newVec.clear();
            oldVec.clear();
            continue;   // should not happen bond has no other bonds
         }

         oldVec.add(oldNat);
         newVec.add(newAtomsAtOldPos[oldNat.GetIdx()]);

         // now set the stereochemistry
         int stereo = oldBd.GetStereo(oldVec, OEBondStereo.CisTrans);
         newBd.SetStereo(newVec, OEBondStereo.CisTrans, stereo);

         newVec.clear();
         oldVec.clear();
      }
      bIt.delete();
      newVec.delete();
      oldVec.delete();

      oechem.OEFindRingAtomsAndBonds(tmpMol);
      oechem.OEAssignAromaticFlags(tmpMol);

      return tmpMol;
   }

   /**
    * Create canonical Smiles adding additional checks to the openeye functionality
    * in order to fix shortcomings:
    *
    * - If mol has stereo atoms invert stereo centers and try atom atom match, if
    *   that works return lexically smaller smiles to circumvent OEBug on eg.
    *   1,4 Dimetyl cyclohexan.
    * - Invert all stero centers non no-chiral atoms and try atom atom match if
    *   that works return lexically smaller smiles.
    * - invert all double bonds if atom atom match yields same compound return
    *   lexically smaller smiles.
    * - generate NCHECK molecule objects with random atom and bond orderings,
    *   compare for smallest smiles.
    */
   public static String myCanonicalSmi(OEGraphMol mol) {
      oechem.OESuppressHydrogens(mol,false,false,true);
      oechem.OEPerceiveChiral(mol);

      int smiFlavor = outFlavor;
      String minSmi = oechem.OECreateSmiString(mol, smiFlavor);

      if( (smiFlavor & OESMILESFlag.ISOMERIC) == 0) return minSmi;

      OEQMol qMol = new OEQMol(mol);
      qMol.BuildExpressions(OEExprOpts.ExactAtoms, OEExprOpts.ExactBonds);
      OESubSearch subSearch = new OESubSearch(qMol);
      qMol.delete();

      // invert all atoms with stereochemistry and check for lowest smiles
      // catches C[C@H]1CCC[C@H]1C
      OEGraphMol tmpMol = new OEGraphMol(mol);
      OEUnaryAtomPred functor = new AtomHasStereoSpecifiedFunctor();
      if( invertAtoms(functor, tmpMol) > 1) {
         if(subSearch.SingleMatch(tmpMol))
         {  // this is a smiles which did not get canonized correctly
            String smi2 = oechem.OECreateSmiString(tmpMol, smiFlavor);
            if(smi2.compareTo(minSmi) < 0 ) {
               minSmi = smi2;
            }
         }
      }
      functor.delete();

      // invert all stereo centers which have stereo but are not chiral
      // catches C[C@H]1CC[C@H](CC1)CC[C@@H](C)Cl
      tmpMol.Clear();
      oechem.OEAddMols(tmpMol, mol);
      functor = new AtomStereoNonChiralFunctor();
      if( invertAtoms(functor, tmpMol) > 1) {
         if(subSearch.SingleMatch(tmpMol))
         {  // this is a smiles which did not get canonized correctly
            String smi2 = oechem.OECreateSmiString(tmpMol, smiFlavor);
            if(smi2.compareTo(minSmi) < 0 ) {
               minSmi = smi2;
            }
         }
      }
      functor.delete();

      // invert all double bonds
      // catches C/C(=N/OC(=O)CCC(=O)O/N=C(\\C)/c1cccs1)/c2cccs2
      tmpMol.Clear();
      oechem.OEAddMols(tmpMol, mol);
      OEUnaryBondPred bFunctor = new BondHasStereoSpecifiedFunctor();
      if( invertBonds(bFunctor, tmpMol) > 0) {
         if(subSearch.SingleMatch(tmpMol))
         {  // this is a smiles which did not get canonized correctly
            String smi2 = oechem.OECreateSmiString(tmpMol, smiFlavor);
            if(smi2.compareTo(minSmi) < 0 ) {
               minSmi = smi2;
            }
         }
      }
      bFunctor.delete();
      subSearch.delete();

      // now generate NCHECK mol objects with random ordering and check for
      // smallest smiles
      for(int i=1; i<=NCHECK; i++) {
         tmpMol = MolReorder.reorder(mol);

         String smi2 = oechem.OECreateSmiString(tmpMol, smiFlavor);
         if( smi2.compareTo(minSmi) < 0 ) {
            minSmi = smi2;
         }
      }
      tmpMol.delete();

      return minSmi;
   }

   /** Invert double bond stereochemistry on bonds in bondFunctor */
   private static int invertBonds(OEUnaryBondPred bondFunctor, OEMolBase mol) {
      int count =0;
      OEAtomBaseVector vec = new OEAtomBaseVector();
      OEBondBaseIter bIt = mol.GetBonds(bondFunctor);
      while( bIt.hasNext() ) {
         OEBondBase bd = bIt.next();
         if(! bd.HasStereoSpecified() ) continue;
         if( bd.IsInRing() ) continue;

         OEAtomBase at1 = bd.GetBgn();
         OEAtomBase at2 = bd.GetEnd();
         vec.add(Atom.findExoNeighbor(at1, at2));
         vec.add(Atom.findExoNeighbor(at2, at1));
         int stereo = bd.GetStereo(vec, OEBondStereo.CisTrans);
         stereo = stereo == OEBondStereo.Cis ? OEBondStereo.Trans : OEBondStereo.Cis;
         bd.SetStereo(vec, OEBondStereo.CisTrans, stereo);
         vec.clear();

         count++;
      }
      vec.delete();
      bIt.delete();

      return count;
   }


   /** Invert atoms identified by atomFunctor */
   private static int invertAtoms(OEUnaryAtomPred atomFunctor, OEMolBase mol) {
      int count =0;
      OEAtomBaseVector atVec = new OEAtomBaseVector();
      OEAtomBaseIter aIt = mol.GetAtoms(atomFunctor);
      while( aIt.hasNext() ) {
         OEAtomBase at = aIt.next();
         invertStereo(at, atVec);
         count++;
      }
      aIt.delete();
      atVec.delete();

      return count;
   }


   private static void invertStereo(OEAtomBase at, OEAtomBaseVector atVec) {
      atVec.clear();
      OEAtomBaseIter neighborsIt=at.GetAtoms();
      while( neighborsIt.hasNext() ) {
         OEAtomBase atNeighbor = neighborsIt.next();
         atVec.add(atNeighbor);
      }
      neighborsIt.delete();

      int stereo = at.GetStereo(atVec, OEAtomStereo.Tetrahedral);
      if(stereo == OEAtomStereo.Left )
         stereo = OEAtomStereo.Right;
      else
         stereo = OEAtomStereo.Left;

      at.SetStereo(atVec, OEAtomStereo.Tetrahedral, stereo);
   }

   public static final void main(String...args)
   {  OEGraphMol mol = new OEGraphMol();
      String smi = "C[C@@H]2[C@@H]([C@H]([C@H]2Cl)C)Cl";
      smi = "C[C@@H]2[C@H]([C@H]([C@@H]2Cl)C)Cl";
      smi = "C[C@@H](O)CCCCCC[C@H](O)C";
      oechem.OEParseSmiles(mol, smi);
      for(int i=0; i<2000; i++) {
         mol = reorder(mol);
         System.err.printf("%4d %s\n",i,myCanonicalSmi(mol));
      }
   }
}
