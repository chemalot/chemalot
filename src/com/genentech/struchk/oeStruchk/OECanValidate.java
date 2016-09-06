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

import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import openeye.oechem.OEAtomBase;
import openeye.oechem.OEAtomBaseIter;
import openeye.oechem.OEAtomBaseVector;
import openeye.oechem.OEAtomStereo;
import openeye.oechem.OEBondBase;
import openeye.oechem.OEBondBaseIter;
import openeye.oechem.OEBondStereo;
import openeye.oechem.OEExprOpts;
import openeye.oechem.OEGraphMol;
import openeye.oechem.OEQMol;
import openeye.oechem.OESMILESFlag;
import openeye.oechem.OESubSearch;
import openeye.oechem.OEUnaryAtomPred;
import openeye.oechem.OEUnaryBondPred;
import openeye.oechem.oechem;
import openeye.oechem.oemolistream;

import com.genentech.oechem.tools.Atom;
import com.genentech.oechem.tools.AtomHasStereoSpecifiedFunctor;
import com.genentech.oechem.tools.AtomIntPropertyFunctor;

/**
 * This class provides a function to randomly reorder atoms in an OEGrapMol.
 *
 * This is to be used to validate the smiles canonicalization routine in oechem.
 * @author albertgo
 *
 */
public class OECanValidate {
   private OECanValidate() {
   }

   /** tag to mark atoms in oechem according to their index */
   public static final int ATOMNumTag = oechem.OEGetTag("STRCHKANum");

   private static final int NCHECK = 30;

   private static final boolean DOWorkarounds = false;

   public static void main(String [] args) throws IOException {
      int outFlavor = OESMILESFlag.AtomStereo | OESMILESFlag.BondStereo
         | OESMILESFlag.Isotopes | OESMILESFlag.Hydrogens;

      System.out.printf("NChecked\tNDiff\tOrgcanSmi\tCanSmi\tIntermediate\n");

      long start = System.currentTimeMillis();

      oemolistream ifs = new oemolistream(args[0]);

      OEGraphMol mol = new OEGraphMol();

      int count = 0;
      Runtime rt = Runtime.getRuntime();
      while ( oechem.OEReadMolecule(ifs , mol ) ) {
         count ++;

         oechem.OESuppressHydrogens(mol,false,false,true);
         oechem.OEPerceiveChiral(mol);
         oechem.OEMDLStereoFromBondStereo(mol);

         if( ! checkChirality(mol) )
         {  System.err.printf("Smiles with invalid stereochenter skiped: %s\n",
               myCanonicalSmi(mol, outFlavor));
            continue;
         }

         checkMol(outFlavor, mol);

         if(count % 1 == 0) System.gc();
         if(count % 100 == 0) System.err.print(".");
         if(count % 4000 == 0)
            System.err.printf( " %d t=%d f=%d u=%d\n",
                count, rt.totalMemory(), rt.freeMemory(), rt.totalMemory() - rt.freeMemory() );

      }
      mol.delete();

      ifs.close();

      System.err.printf("%d compounds in %d sec\n", count,  (System.currentTimeMillis()-start)/3600);
   }

   /**
    * Generate NCHECK random orders of atoms and compare canonical smiles to
    * canonical smiles of input.
    */
   private static int checkMol(int outFlavor, OEGraphMol mol) {

      int countDiffs = 0;
      // create canonical smiles as basis.
//      String orgCanSmi = oechem.OECreateSmiString(mol, outFlavor|OESMILESFlag.Canonical);
      // myCanonicalSmi already checks for some obvious cases
      String orgCanSmi = myCanonicalSmi(mol, outFlavor);

      // create NCHECK random smiles and check cannonization
      Set<String> checkedSmis = new HashSet<String>();
      Set<String> printedSmis = new HashSet<String>();
      for(int i=1; i<=NCHECK; i++) {
         OEGraphMol outMol = MolReorder.reorder(mol);

         oechem.OESuppressHydrogens(outMol,false,false,true);

         String smi = oechem.OECreateSmiString(outMol, outFlavor);
         if(checkedSmis.contains(smi)) {
            outMol.delete();
            continue;
         }
         checkedSmis.add(smi);

         outMol.Clear();
         oechem.OEParseSmiles(outMol, smi);
         oechem.OESuppressHydrogens(outMol,false,false,true);
         String canSmi = myCanonicalSmi(outMol, outFlavor);
         outMol.delete();

         if(! orgCanSmi.equals(canSmi) && ! printedSmis.contains(canSmi)) {
            countDiffs++;
            System.out.printf("\n%d %d %s\t%s\t%s\n", i, countDiffs, orgCanSmi, canSmi, smi);
            printedSmis.add(canSmi);
            System.out.flush();
         }
      }
      return countDiffs;
   }

   /**
    * If mol has stereo atoms invert stereo centers and try atom atom match, if
    * that works return lexically smaller smiles to circumvent OEBug on eg.
    * 1,4 Dimetyl cyclohexan.
    */
   private static String myCanonicalSmi(OEGraphMol mol, int smiFlavor) {
      smiFlavor = smiFlavor|OESMILESFlag.Canonical;
      String smi1 = oechem.OECreateSmiString(mol, smiFlavor);
      if( (smiFlavor & OESMILESFlag.ISOMERIC) == 0) return smi1;

      if( DOWorkarounds )
      {  oechem.OEPerceiveChiral(mol);

         // invert all atoms with stereochemistry and check for lowest smiles
         // catches C[C@H]1CCC[C@H]1C
         OEUnaryAtomPred functor = new AtomHasStereoSpecifiedFunctor();
         smi1 = invertAtomAndCheck(mol, functor, smi1, smiFlavor);
         functor.delete();

         // catches C[C@H]1CC[C@H](CC1)CC[C@@H](C)Cl
         functor = new AtomStereoNonChiralFunctor();
         smi1 = invertAtomAndCheck(mol, functor, smi1, smiFlavor);

         // catches C/C(=N/OC(=O)CCC(=O)O/N=C(\\C)/c1cccs1)/c2cccs2
         OEUnaryBondPred bFunctor = new BondHasStereoSpecifiedFunctor();
         smi1 = invertBondAndCheck(mol, bFunctor, smi1, smiFlavor);
         bFunctor.delete();
      }

      return smi1;
   }

   /**
    * Invert stereo centers returned by atomFunctor and check if molecule is identical
    * after inversion, if identical return lexically smaller smiles.
    * @param mol molecule object on which to perform check
    * @param atomFunctor functor returning atoms in mol which are to be inverted
    * @param orgSmiles canonical smiles of mol
    * @param smiFlavor canonicalization flavor to be used
    * @return either orgSmiles or smiles after inversion if in yeilds same molecule with smaller smiles
    */
   private static String invertAtomAndCheck(OEGraphMol mol,
         OEUnaryAtomPred atomFunctor, String orgSmiles, int smiFlavor) {
      OEQMol invertMol = new OEQMol(mol);

      // invert all atoms with hasStereoSpecified and check for identity
      // use smaller smiles if identical

      // one stereo or less stereo centers can not have bug
      if( invertAtoms(atomFunctor, invertMol) <= 1) {
         invertMol.delete();
         return orgSmiles;
      }

      invertMol.BuildExpressions(OEExprOpts.ExactAtoms, OEExprOpts.ExactBonds);
      OESubSearch subSearch = new OESubSearch(invertMol);
      // no match => mirror image is enantiomer
      if(! subSearch.SingleMatch(mol)) {
         subSearch.delete();
         invertMol.delete();
         return orgSmiles;
      }
      subSearch.delete();

      // this is a smiles which did not get canonized correctly
      String smi2 = oechem.OECreateSmiString(invertMol, smiFlavor);
      invertMol.delete();

      if(smi2.compareTo(orgSmiles) < 0) orgSmiles = smi2;

      return orgSmiles;
   }

   /**
    * Invert stereo double bonds returned by bondFunctor and check if molecule is identical
    * after inversion, if identical return lexically smaller smiles.
    * @param mol molecule object on which to perform check
    * @param bondFunctor functor returning bonds in mol which are to be inverted
    * @param orgSmiles canonical smiles of mol
    * @param smiFlavor canonicalization flavor to be used
    * @return either orgSmiles or smiles after inversion if in yeilds same molecule with smaller smiles
    */
   private static String invertBondAndCheck(OEGraphMol mol,
         OEUnaryBondPred bondFunctor, String orgSmiles, int smiFlavor) {
      OEQMol invertMol = new OEQMol(mol);

      // invert all atoms with hasStereoSpecified and check for identity
      // use smaller smiles if identical

      // one stereo or less stereo centers can not have bug
      if( invertBonds(bondFunctor, invertMol) < 1) {
         invertMol.delete();

         return orgSmiles;
      }

      invertMol.BuildExpressions(OEExprOpts.ExactAtoms, OEExprOpts.ExactBonds);
      OESubSearch subSearch = new OESubSearch(invertMol);
      // no match => mirror image is enantiomer
      if(! subSearch.SingleMatch(mol)) {
         invertMol.delete();
         subSearch.delete();
         return orgSmiles;
      }
      subSearch.delete();

      // this is a smiles which did not get canonized correctly
      String smi2 = oechem.OECreateSmiString(invertMol, smiFlavor);

      if(smi2.compareTo(orgSmiles) < 0) orgSmiles = smi2;

      return orgSmiles;
   }

   private static int invertBonds(OEUnaryBondPred bondFunctor, OEQMol mol) {
      int count =0;
      OEBondBaseIter bIt = mol.GetBonds(bondFunctor);
      while( bIt.hasNext() ) {
         OEBondBase bd = bIt.next();
         if(! bd.HasStereoSpecified() ) continue;
         if( bd.IsInRing() ) continue;

         OEAtomBase at1 = bd.GetBgn();
         OEAtomBase at2 = bd.GetEnd();
         OEAtomBaseVector vec = new OEAtomBaseVector();
         vec.add(Atom.findExoNeighbor(at1, at2));
         vec.add(Atom.findExoNeighbor(at2, at1));
         int stereo = bd.GetStereo(vec, OEBondStereo.CisTrans);
         stereo = stereo == OEBondStereo.Cis ? OEBondStereo.Trans : OEBondStereo.Cis;
         bd.SetStereo(vec, OEBondStereo.CisTrans, stereo);
         vec.delete();

         count++;
      }
      bIt.delete();

      return count;
   }

   private static int invertAtoms(OEUnaryAtomPred atomFunctor, OEQMol mol) {
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


   private static void invertStereo(OEAtomBase at) {
      OEAtomBaseVector atVec = new OEAtomBaseVector();
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
      atVec.delete();
   }

   static final String getAtBd(OEGraphMol inMol) {
      StringBuilder sb = new StringBuilder();

      OEAtomBaseIter nAtIt=inMol.GetAtoms();
      while( nAtIt.hasNext() ) {
         OEAtomBase at= nAtIt.next();
         sb.append(oechem.OEGetAtomicSymbol(at.GetAtomicNum()));
      }
      nAtIt.delete();

      OEBondBaseIter bIt = inMol.GetBonds();
      while( bIt.hasNext() ) {
         OEBondBase bd = bIt.next();
         sb.append(bd.GetOrder());
      }
      bIt.delete();

      return sb.toString();
   }


   /** return true if atom has atoms which have sterochemsitry (wedge bonds)
    * but are not chiral and do not need wedge bonds.
    *
    * This will return true for C[C@H](C)C but not for C[C@H]1C[C@H]1C.
    *
    */
   private static boolean checkChirality(OEGraphMol mol) {
      // tag atoms so that we can distinguish them afterwards
      int aNum = 0;
      OEAtomBaseIter aIt = mol.GetAtoms();
      while( aIt.hasNext() )
         aIt.next().SetIntData(ATOMNumTag, aNum++);
      aIt.delete();

      OEQMol qMol = new OEQMol(mol);
      qMol.BuildExpressions(OEExprOpts.ExactAtoms, OEExprOpts.ExactBonds);
      OESubSearch subSearch = new OESubSearch(qMol);
      qMol.delete();

      aIt = mol.GetAtoms();
      while( aIt.hasNext() ) {
         OEAtomBase at = aIt.next();

         if(at.HasStereoSpecified() && ! at.IsChiral() ) {
            if(checkAChiralStereoAtom(subSearch, mol, at)) {
               aIt.delete();
               subSearch.delete();
               return false;
            }
         }
      }
      aIt.delete();
      subSearch.delete();

      return true;
   }

   private static boolean checkAChiralStereoAtom(OESubSearch subSearch, OEGraphMol in, OEAtomBase at) {

      OEGraphMol tmpMol = new OEGraphMol(in);

      // get atom on copy
      AtomIntPropertyFunctor functor =
         new AtomIntPropertyFunctor(OEStruchk.ATOMNumTag,
                                          at.GetIntData(OEStruchk.ATOMNumTag));
      OEAtomBase tmpAt = tmpMol.GetAtom(functor);
      functor.delete();

      OEAtomBaseVector v = new OEAtomBaseVector();
      OEAtomBaseIter nbriter=tmpAt.GetAtoms();
      while( nbriter.hasNext() )
         v.add(nbriter.next());
      nbriter.delete();

      // get current stereo;
      int stereo = tmpAt.GetStereo(v, OEAtomStereo.Tetrahedral);
      // invert
      stereo = (stereo == OEAtomStereo.Left ? OEAtomStereo.Right : OEAtomStereo.Left);
      tmpAt.SetStereo(v, OEAtomStereo.Tetrahedral, stereo);
      v.delete();

      boolean isACiral = false;
      if(subSearch.SingleMatch(tmpMol)) isACiral = true;

      tmpMol.delete();
      return isACiral;
   }
}


