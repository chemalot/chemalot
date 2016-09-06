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

import java.util.*;

import openeye.oechem.*;

import org.jdom.Element;

import com.aestel.chemistry.depict.DepictHelper;
import com.aestel.chemistry.depict.ImageType;
import com.aestel.utility.MessageList;
import com.genentech.oechem.tools.Atom;
import com.genentech.oechem.tools.AtomIntPropertyFunctor;
import com.genentech.struchk.oeStruchk.OEStruchk.StructureFlag;

public class FlagNonChiralStereoCenters extends AbstractStructureCheck
{  private final OEGraphMol copyAllStereoSpecified = new OEGraphMol();
   private final OEGraphMol tmpMol2 = new OEGraphMol();
   private final OEAtomBaseVector atVec = new OEAtomBaseVector();

   private final StrainedCycleStereoGenerator[] strainedCycle;

   @SuppressWarnings("unchecked")
   FlagNonChiralStereoCenters(Element elem)
   {  super(elem );

      ArrayList<StrainedCycleStereoGenerator> cycleList = new ArrayList<StrainedCycleStereoGenerator>();
      for(Element cycleEle : ((List<Element>)elem.getChildren("strainedBicycle")))
      {  String smiles = cycleEle.getTextTrim();
         cycleList.add(new StrainedCycleStereoGenerator(smiles));
      }
      strainedCycle = cycleList.toArray(new StrainedCycleStereoGenerator[cycleList.size()]);
   }



   @Override
   public boolean checkStructure(OEGraphMol in, StructureFlag inStereo,
            MessageList msgs)
   {  // prevent duplicate execution for same molecule
      // because other checks might be resetting this flag this
      // will be executed twice in the checking procedure for each molecule
      // one initial time to make sure that the invalid atoms are tagged
      // and once in the stereoflag assignment/checking
      if( in.GetBoolData(OEStruchk.NONChiralStereoAssignedTag))
         return true;

      in.SetBoolData(OEStruchk.NONChiralStereoAssignedTag, true);

      flagNonchiralStereoAtoms(in);

      return true;
   }


   @Override
   public String getDescriptionHTML()
   {  StringBuilder sb = new StringBuilder(2000);
      sb.append(getDescription()).append("<br/>List of Strained ring systems:\n<table><tr>");

      for( int row=0; row < (strainedCycle.length-1)/4+1; row ++)
      {  for( int i=0; i<4 && row*4+i < strainedCycle.length; i++)
         {  sb.append("<td>")
              .append(DepictHelper.AESTELDepictHelper.getPublicSmilesImageElement(
                       80,100,ImageType.PNG,strainedCycle[row*4+i].getSmiles()))
              .append("</td>");
         }
         sb.append("</tr><tr>\n");

         for( int i=0; i<4 && row*4+i < strainedCycle.length; i++)
         {  sb.append("<td>")
              .append(row*4+i).append("<br/><font size='1'>")
              .append(AbstractStructureCheck.encodeHTML(strainedCycle[row*4+i].getSmiles()))
              .append("</font></td>");
         }
         sb.append("</tr><tr>\n");
      }
      sb.setLength(sb.length()-5); // remove last tr
      sb.append("</table>\n");

      return sb.toString();
   }



   @Override
   public void delete()
   {  atVec.delete();
      copyAllStereoSpecified.delete();
      tmpMol2.delete();
      for( StrainedCycleStereoGenerator cy : strainedCycle)
         cy.delete();
   }

   /**
    * Flag stereogenic sp3 atoms {@link #NONChiralStereoAtomTag} if they
    * they do contain stereo information.
    *
    * Examples:
    *    C[C@H]1CC[C@@H](C)CC1
    *    C[C@H](F)[C@H](C)[C@H](F)C
    *
    * It will also flag some atoms in highly symmetrical cage system such as:
    * C[C@]23C[C@](C[C@H]4C3)(F)C[C@@H](C4)C2.
    * as non chiral.
    *
    * The basic principal involved here is that if a stereo center is inverted
    * and it yields to a graph which can not be overlaid with the input the
    * center is stereogenic. However there are complications due to strained
    * ring systems, and interactions of multiple centers which are explained
    * below in the code.
    */
   private void flagNonchiralStereoAtoms(OEGraphMol in) {

      oechem.OEPerceiveChiral(in);
      flagChiralNitrogens(in);
      makeCopyWithStrainedStereoDefined(in);

      Set<OEAtomBase> candidateSet = new HashSet<OEAtomBase>();

      int[] nSymClassByAtIdx = new int[in.GetMaxAtomIdx()+1];
      int[] nRingBondsByAtomIdx = new int[in.GetMaxAtomIdx()+1];

      // find candidates which have 3 symmetry classes
      int nChiral = 0;
      int nChiralInRing = 0;
      int nSymSpiro = 0;
      OEAtomBaseIter aIt = in.GetAtoms();
      while(aIt.hasNext()) {
         OEAtomBase at = aIt.next();
         if( at.IsAromatic() ) continue; // Cc1cs(=O)cc1

         if( at.IsChiral() || at.GetBoolData(OEStruchk.ISChiralNotRecognized) ) {
            nChiral++;
            if( at.IsInRing() )
               nChiralInRing++;
            continue;
         }

         if( at.IsCarbon() ) {   // carbon
            if( at.GetExplicitDegree() < 3 ) continue;   // CH2 or CH3
            if( at.GetDegree() != 4 ) continue;         // not thetraedral carbon

         }else if( at.IsSulfur() ) { // Sulfur C[S@](=O)CC
            // we are only supporting stereo on *-S(=X)-* and S(=O)(=N)??? should we support -S(=X)-H???
            if( at.GetExplicitDegree() < 3 ) continue;

            // this is not chiral C[S+](C[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@@H]([C@@H]1O)O)CC[C@H](N)C(O)=O
            if( at.GetFormalCharge() != 0 ) continue;

         }else if( at.IsNitrogen() ) { // OC1C[N+]2(C1)CC[N+]3(CC2)CC(O)C3
            if( at.GetExplicitDegree() != 4 ) continue;
            if( at.GetFormalCharge() != 1 ) continue;
            if( at.GetTotalHCount() != 0 ) continue;

         } else {
            continue; // we expect non chiral stereogenic center on C and S only
         }


         int nRingBonds = getNRingBonds(at);
         nRingBondsByAtomIdx[at.GetIdx()] = nRingBonds;

         int nSymClass  = getNSymClass(at);
         nSymClassByAtIdx[at.GetIdx()] = nSymClass;

         if(    nSymClass <= 2
             && nRingBonds == 4
             && ! at.GetBoolData(OEStruchk.ISStrainedBridgeHead) )
         {  // spiro carbon needs to keep stereo in: C[C@H]1CC[C@]2(CC1)CC[C@@H](CC2)C
            // mark for check by inversion
            candidateSet.add(at);
         if( nSymClass == 2 ) nSymSpiro++; // eg. N1C[C@]12CN2

            continue;
         }


         if( nSymClass == 4 ) // Bug in OEChem isChiral: F[C@@H]([3H])C
         {  at.SetChiral(true);
            at.SetBoolData(OEStruchk.ISChiralNotRecognized, true);
         }

         if( nSymClass != 3 ) continue;
         candidateSet.add(at);
      }
      aIt.delete();


      /////////////////////////////////////////////////////////////////////////
      // now we have found all potential candidates for nonchiral stereogenic centers
      // validate results

      // TODO this is too stringent and will caus ether center in N1C[C@]12CN2 to be missed
      if( nChiralInRing == 0 )
      {  if( candidateSet.size() == 0 ) return;
         if( candidateSet.size() == 1              // not CC1CC(C)C1
             && nSymSpiro == 0                     // not N1C[C@]12CN2
             && nChiral < 2                        // not C[C@H](C[C@@H](C)F)C[C@H](C)F
             && getNStereoBonds(in) < 2) return;   // C/C=C/[C@H](C)/C=C\C
      }

      int[] ringSysByAtom = new int[in.GetMaxAtomIdx()];
      int nRingSys = oechem.OEDetermineRingSystems(in, ringSysByAtom);

      Set<Integer> candidateNumSet = new HashSet<Integer>(candidateSet.size());

      for( OEAtomBase at : candidateSet)
         candidateNumSet.add(at.GetIntData(OEStruchk.ATOMNumTag));
      assignStereoToCandidatesInCopy(candidateNumSet);


      // now process all candidates by trying to invert stereo
      // and seeing if they still match the input
      if( candidateSet.size() > 0 )
      {
         // create Substructure search where all candidates have stereo specified
         OESubSearch subSearchAllSpecified = buildSubSearch(copyAllStereoSpecified);

         // try to invert non-bridgehead candidate atoms one by one
         // if that results in a non-matchable structure this is a real NONChiralStereoAtom
         for( OEAtomBase at : candidateSet)
         {  if(    ! at.GetBoolData(OEStruchk.ISStrainedBridgeHead)
                && isAChiralStereo(subSearchAllSpecified, at) )
            {  at.SetBoolData(OEStruchk.NONChiralStereoAtomTag, true);
               // N's in S[C@@H](O1)O[N@+]21OO[N@@+]3(O[C@@H](S)O3)OO2
               // where flagged as not having stereo
               at.DeleteData(OEStruchk.STEREOClearTag);
            }
         }
         subSearchAllSpecified.delete();
      }

      // now if we found at least one NONChiralStereoAtom
      // we might have a situation like: C1[C@@H](C[C@H](C[C@@H]1F)F)F
      // in this case inversion of two of the stereo centers will have resulted
      // in a match so the two atoms will not have been flagged as stereogenic.
      // Lets go back and in all ring system flag all candidates if there is at
      // least one detected NONChiralStereoAtom
      for(int ring=1; ring <= nRingSys; ring++)
      {  List<OEAtomBase> ringCandidates = getRingCandidates(in, ring, ringSysByAtom, candidateSet);
         /* eg 13C in O[13CH]1C(C2)CC3CC1CC2C */
         List<OEAtomBase> ringNonStrained = new ArrayList<OEAtomBase>();

         /* eg all other bridgehead in O[13CH]1C(C2)CC3CC1CC2C */
         List<OEAtomBase> ringStrainedHeads = new ArrayList<OEAtomBase>();

         for(int i=0; i<ringCandidates.size(); i++)
         {  OEAtomBase at = ringCandidates.get(i);

            if( at.GetBoolData(OEStruchk.ISStrainedBridgeHead) )
               ringStrainedHeads.add(at);
            else if( at.GetBoolData(OEStruchk.NONChiralStereoAtomTag) )
               ringNonStrained.add(at);
         }

         if( ringNonStrained.size() == 1 )
         {  if( ringStrainedHeads.size() > 0)
            {  // could have no stereo as O[13CH]1C(C2)CC3CC1CC2C3
               // could have stereo as    O[13C@H]1[C@@H](C2)C[C@@H]3C[C@H]1C[C@@]2(C)C3
               // check to see if stereo stays even when all bridgehead atoms have been
               //   removed of parity
               OEAtomBase at = ringNonStrained.get(0);
               if( ! invertBridgeAndCheckStereogenic(in, at, ringStrainedHeads))
               {  at.SetBoolData(OEStruchk.STEREOClearTag, true);
                  OEAtomBase tmpAt = getAtomByNum(copyAllStereoSpecified, at);
                  Atom.removeChiralInfo(tmpAt);
                  ringNonStrained.clear();
               }
            }
         }

         int nCleared = 0;
         if( ringNonStrained.size() > 0 )
         {  // if there is one non chiral sterogenic center in a ring system such as
            // O[13C@H]1[C@@H](C2)C[C@@H]3C[C@H]1C[C@@]2(C)C3
            // all strained bridge heads also are stereogenic
            for( OEAtomBase at : ringCandidates )
            {  if( at.HasData( OEStruchk.NONChiralStereoAtomTag) ) continue; // already assigned

               if( symmetricSubstituentsHaveStereo(in, at, candidateNumSet ) )
               {  // O[13C@H]1[C@@H](C2)C[C@@H]3C[C@H]1C[C@@]2(C)C3
                  // C1[C@@H](C[C@@H](C[C@H]1F)F)F
                  at.SetBoolData(OEStruchk.NONChiralStereoAtomTag, true);

               }else
               {  // the dimethyl in: CC1([C@]2(C(CC=C2)(C)C)C=CC1)C
                  at.SetBoolData(OEStruchk.STEREOClearTag, true);
                  OEAtomBase tmpAt = getAtomByNum(copyAllStereoSpecified, at);
                  Atom.removeChiralInfo(tmpAt);
                  nCleared++;
               }
            }
         }

         if( ringNonStrained.size() - nCleared == 0 ) // no exeocyclic stereogenic left
         {  // -> all bridge heads which are nonChiral are alos not stereogenic
            for( OEAtomBase at : ringStrainedHeads )
            {  OEAtomBase tmpAt = getAtomByNum(copyAllStereoSpecified, at);
               Atom.removeChiralInfo(tmpAt);
            }
         }

         List<OEAtomBase> ringChiral = getRingChiral(in, ring, ringSysByAtom);
         invertCagedBridgeHeadsAndCheckChiral(in, ringChiral);
      }
   }


   /**
    * Check if OC bond is stereogenic in O[13C@H]1[C@@H](C2)C[C@@H]3C[C@H]1C[C@@]2(C)C3.
    *
    * IF after removing parity from all strained bridgehead atoms inversion
    * of the exoAtom not be matched onto itself the atom is stereogenic.
    * @param orgExoAt
    * @param ringStrainedHeads straiend bridgeheads in same rings system
    * @return
    */
   private boolean invertBridgeAndCheckStereogenic(OEGraphMol in,
                             OEAtomBase orgExoAt, List<OEAtomBase> ringStrainedHeads)
   {  tmpMol2.Clear();
      oechem.OEAddMols(tmpMol2, copyAllStereoSpecified);

      for(OEAtomBase orgAt : ringStrainedHeads)
      {  if( orgAt.IsChiral() ) continue;

         OEAtomBase tmpAt = getAtomByNum(tmpMol2, orgAt);
         Atom.removeChiralInfo(tmpAt);
      }

      OEAtomBase exoAt = getAtomByNum(tmpMol2, orgExoAt);
      OESubSearch ss = buildSubSearch(tmpMol2);
      invertStereo(exoAt);

      boolean isStereogenic = true;
      if( ss.SingleMatch(tmpMol2) ) isStereogenic = false;
      ss.delete();

      return isStereogenic;
   }


   /** get list of chiral atoms in ring */
   private static List<OEAtomBase> getRingChiral(OEGraphMol in, int ring, int[] ringSysByAtom)
   {  ArrayList<OEAtomBase> chList = new ArrayList<OEAtomBase>();

      OEIsChiralAtom isCA = new OEIsChiralAtom();
      OEAtomBaseIter atIt = in.GetAtoms(isCA);
      while(atIt.hasNext())
      {  OEAtomBase at= atIt.next();
         if( ringSysByAtom[at.GetIdx()] == ring)
            chList.add(at);
      }
      atIt.delete();

      return chList;
   }


   private static List<OEAtomBase> getRingCandidates(OEGraphMol in, int ring,
                                              int[] ringSysByAtom, Set<OEAtomBase> candidateList)
   {  ArrayList<OEAtomBase> rcList = new ArrayList<OEAtomBase>();
      for( OEAtomBase a : candidateList )
         if( ringSysByAtom[a.GetIdx()] == ring )
            rcList.add(a);

      return rcList;
   }


   /** get the number of double bonds with cis trans stereochemistry */
   private static int getNStereoBonds(OEGraphMol in)
   {  int nChiralBonds = 0;
      OEBondBaseIter bIt = in.GetBonds();
      while(bIt.hasNext())
      {  OEBondBase bd = bIt.next();
         if( bd.GetOrder() == 2 && bd.IsChiral() && bd.HasStereoSpecified(OEBondStereo.CisTrans))
            nChiralBonds++;
      }
      bIt.delete();

      return nChiralBonds;
   }


   /** create a copy of in into copyAllStereoSpecified which has all bridgehead
    * stereo atom defined with acceptable but random stereochemistry.
    *
    * Also, flag strained bridge heads.
    *
    */
   private void makeCopyWithStrainedStereoDefined(OEMolBase in)
   {  copyAllStereoSpecified.Clear();
      oechem.OEAddMols(copyAllStereoSpecified, in);
      oechem.OESuppressHydrogens(copyAllStereoSpecified,false,false,true);

      // define stereo on caged ring systems
      for( StrainedCycleStereoGenerator cy : strainedCycle )
         cy.assignStereo(copyAllStereoSpecified, in);
   }



   /** for all candidates assign a random stereo chemistry in the
    * copyAllStereoSpecified of the current molecule.
    *
    * This is used then later when inverting single centers to check if the center is
    * a real stereogenic center.
    */
   private void assignStereoToCandidatesInCopy(Set<Integer> candidateSet )
   {  // invent stereo for any other candidates or chiral atoms
      OEAtomBaseIter atIt = copyAllStereoSpecified.GetAtoms();
      while(atIt.hasNext() )
      {  OEAtomBase at = atIt.next();
         int atNum = at.GetIntData(OEStruchk.ATOMNumTag);

         if(    candidateSet.contains(atNum)
               || (at.IsChiral()
                   && ! at.GetBoolData(OEStruchk.STEREOClearTag) // isChiral has not been found invalid
                   && (at.IsCarbon()||at.IsSulfur()||at.IsPhosphorus()||at.IsNitrogen())))
//               || (at.IsChiral() && (at.IsCarbon()||at.IsSulfur()||at.IsPhosphorus()))
//               || (at.IsNitrogen() && at.GetBoolData(OEStruchk.ISChiralNotRecognized)) )
         {  if( at.HasStereoSpecified() ) continue; // no need to define

            getStereoNeighbors(at, atVec);

            int stereo = OEAtomStereo.Right;
            if( ! at.SetStereo(atVec, OEAtomStereo.Tetrahedral,stereo))
           {  stereo = OEAtomStereo.Left;
               at.SetStereo(atVec, OEAtomStereo.Tetrahedral,stereo);
            }
         }
      }
      atIt.delete();
      oechem.OEPerceiveChiral(copyAllStereoSpecified); // needed???
   }


   /** get the neighbor atoms of at into atVec2.
    *
    * If possible make sure the first atom is not in a ring.
    */
   private static void getStereoNeighbors(OEAtomBase at, OEAtomBaseVector atVec2)
   {  atVec2.clear();
      OEAtomBaseIter atIt=at.GetAtoms();
      while( atIt.hasNext() )
      {  OEAtomBase at2 = atIt.next();
         atVec2.add(at2);
      }
      atIt.delete();
   }


   /** flag all bridgehead nitrogens with 3 symmetry classes.
    * @param nRingBondsByAtIdx */
   private static void flagChiralNitrogens(OEGraphMol in)
   {  OEIsNitrogen npred = new OEIsNitrogen();
      OEAtomBaseIter ait = in.GetAtoms(npred);
      while( ait.hasNext() )
      {  OEAtomBase at = ait.next();
         if(! isChiralNitrogen(at) )
         {  at.SetChiral(false);
            at.SetBoolData(OEStruchk.STEREOClearTag, true);
         }
      }
      ait.delete();
      npred.delete();
   }

   private static boolean isChiralNitrogen(OEAtomBase at)
   {  assert at.IsNitrogen();

      if( !at.IsChiral() ) return false;
      if( at.GetValence() > at.GetDegree() ) return false; // has double bond
      if( at.IsAromatic() ) return false;
      if( at.GetExplicitDegree() < 3 ) return false;
      if( at.GetExplicitDegree() == 4 && at.GetTotalHCount() == 0 ) return true;

      if( getNRingBonds(at) < 3 ) return false;

      if( ! at.GetBoolData(OEStruchk.ISStrainedBridgeHead)) return false;

      int nBridgeHeadNeighbors = 0;
      OEAtomBaseIter atIt = at.GetAtoms();
      while(atIt.hasNext())
      {  OEAtomBase n = atIt.next();
         if( n.GetBoolData(OEStruchk.ISStrainedBridgeHead) || n.IsAromatic())
            nBridgeHeadNeighbors++;
      }
      atIt.delete();

      // if there is one or more bridgehead neighbors assume this is invertible
      if( nBridgeHeadNeighbors >= 1) return false;

      return true;
   }


   static boolean canBeChiral(OEAtomBase at)
   {  if( at.IsCarbon()||at.IsSulfur()||at.IsPhosphorus()||at.IsNitrogen())
         return true;
      return false;
   }

   /**
    * Invert atoms stereo and see if it is still matched by a SSS with the original stereo.
    *
    * Assume all candidates have been set to R stereo.
    *
    * @return true if atom is really a achiral stereogenic center.
    */
   private boolean isAChiralStereo( OESubSearch subSearchAllSpecified, OEAtomBase srcAtom ) {
      tmpMol2.Clear();
      oechem.OEAddMols(tmpMol2, copyAllStereoSpecified);

      // get atom on copy
      OEAtomBase tmpAt = getAtomByNum(tmpMol2, srcAtom);

      invertStereo(tmpAt);

      if( subSearchAllSpecified.SingleMatch(tmpMol2) ) return false;

      return true;
   }

   /**
    * Verify that atom with 4 substituents in 3 symentry classes has stereo on the two
    * substituents with the identical symmetry class.
    *
    * If not thN this is a case like: CC1([C@]2(C(CC=C2)(C)C)C=CC1)C and the atom
    * is not stereogenic.
    */
   private static boolean symmetricSubstituentsHaveStereo(OEMolBase mol, OEAtomBase at, Set<Integer> candidateNumSet)
   {  OEAtomBase[] neighbors = new OEAtomBase[4];
      int nNeighbors = 0;

      OEAtomBaseIter atIt = at.GetAtoms();
      while(atIt.hasNext())
         neighbors[nNeighbors++] = atIt.next();
      atIt.delete();

      assert nNeighbors >= 3;

      for(int i=0; i < nNeighbors; i++)
      {  for(int j=i+1; j < nNeighbors; j++)
         {  if( neighbors[i].GetSymmetryClass() == neighbors[j].GetSymmetryClass() )
            {  // i and j are the two symmetry equivalent neighbors
               // do depth first search establishing if there are stereocenters along the
               // bat from at to neighbor[i].
               // if there are we have a case like C1[C@@H](C[C@@H](C[C@H]1F)F)F
               // if there are not we have a case like the dimethyl in: CC1([C@]2(C(CC=C2)(C)C)C=CC1)C
               boolean[] visited = new boolean[mol.GetMaxAtomIdx()+1];
               Arrays.fill(visited, false);
               visited[at.GetIdx()] = true;
               if( findStereoInPath(visited, neighbors[i], candidateNumSet) )
                  return true;
               return false;
            }
         }
   }
      assert false: "we should not be here as the assumption is that this atom has "
                   +"3 symmetry classes and two substitutions which are the same.";
      return true;
   }


   /** Recursively look for atoms which have either stereo or are in candidateNumSet
    *  starting from at.
    */
   private static boolean findStereoInPath(boolean[] visited, OEAtomBase at,
            Set<Integer> candidateNumSet)
   {  if( at.IsChiral() ) return true;
      if( candidateNumSet.contains(at.GetIntData(OEStruchk.ATOMNumTag))) return true;

      visited[at.GetIdx()] = true;

      OEAtomBaseIter atIt = at.GetAtoms();
      while( atIt.hasNext() )
      {  at = atIt.next();
         if( ! visited[at.GetIdx()] )
         {  if( findStereoInPath(visited, at, candidateNumSet) )
            {  atIt.delete();
               return true;
            }
         }
      }
      atIt.delete();
      return false;
   }


   private static int getNRingBonds(OEAtomBase at)
   {  int nRingBonds = 0;
      OEBondBaseIter bIt = at.GetBonds();

      while(bIt.hasNext()) {
         OEBondBase bd = bIt.next();
         if( bd.IsInRing()) nRingBonds++;
      }
      bIt.delete();

      return nRingBonds;
   }


   /** return the number of distinct symmetry classes of the atoms neighbors */
   private static int getNSymClass(OEAtomBase at)
   {  List<AtomBond> atomList = new ArrayList<AtomBond>(4);
      OEBondBaseIter bIt = at.GetBonds();

      // get list of neighbors
      while(bIt.hasNext())
      {  OEBondBase bd = bIt.next();

         AtomBond ab = new AtomBond(bd, bd.GetNbr(at));
         atomList.add(ab);
      }
      bIt.delete();

      // sort by symmetry class
      Collections.sort(atomList, new Comparator<AtomBond>() {
         @Override
         public int compare(AtomBond o1, AtomBond o2) {
            return o2.otherAtom.GetSymmetryClass() - o1.otherAtom.GetSymmetryClass();
         }
      });

      // count symmetry classes
      int nSymClass = 0;
      int lastSymClass = -1;
      for(AtomBond ab: atomList) {
         if(ab.otherAtom.GetSymmetryClass() != lastSymClass)
         {  nSymClass++;
            lastSymClass = ab.otherAtom.GetSymmetryClass();
         }
      }

      if( at.GetImplicitHCount() > 0 )
      {  if( at.GetExplicitHCount() == 0 )
         {  nSymClass++;

         } else
         {  // has both implicit and explicit H, check to see if the isotopes are the same
            OEIsHydrogen isH = new OEIsHydrogen();
            OEAtomBaseIter aIt = at.GetAtoms(isH);
            boolean hasNoIsotopeExplicitH = false;
            while( aIt.hasNext() )
               if( aIt.next().GetIsotope() == 0 )
                  hasNoIsotopeExplicitH = true;
            aIt.delete();
            isH.delete();

            if( ! hasNoIsotopeExplicitH ) nSymClass++;

         }
      }

      if( at.IsSulfur() && at.GetDegree() ==3 && at.GetValence() == 4)
         nSymClass++;  // count lone pair in CS(C)=O

      if( at.IsNitrogen() && at.GetBoolData(OEStruchk.ISStrainedBridgeHead)
          && at.GetDegree() == 3)
         nSymClass++;  // count lone pair in C[C@]12C[N@@](C[C@@]3(C)[C@H]2[O])[C@](F)([H])[N@@](C3)C1

      return nSymClass;
   }


   private void invertStereo(OEAtomBase tmpAt)
   {  getStereoNeighbors(tmpAt, atVec);

      // get current stereo;
      int stereo = tmpAt.GetStereo(atVec, OEAtomStereo.Tetrahedral);
      // invert
      stereo = stereo == OEAtomStereo.Left ? OEAtomStereo.Right : OEAtomStereo.Left;
      tmpAt.SetStereo(atVec, OEAtomStereo.Tetrahedral, stereo);
   }


   /** check for highly symmetrical ring systems where atoms with 4 symmetry classes
    * are not chiral eg. CC12CC(CC(C2)(O)C3)CC3C1 C1(NN2)OOC2CC1
    * @param ringStrainedHeads
    * @param ringNonStrained
    * @param subSearchAllSpecified
    */
   private void invertCagedBridgeHeadsAndCheckChiral(OEMolBase mol, List<OEAtomBase> ringChiral)
   {  if( ringChiral.size() == 0 ) return;

      tmpMol2.Clear();
      oechem.OEAddMols(tmpMol2, copyAllStereoSpecified);

      OESubSearch ss = buildSubSearch(tmpMol2);
      int nInverted = 0;
      for(OEAtomBase orgAt : ringChiral)
      {  if( ! orgAt.GetBoolData(OEStruchk.ISStrainedBridgeHead)) continue;

         OEAtomBase tmpAt = getAtomByNum(tmpMol2, orgAt);
         invertStereo(tmpAt);
         nInverted++;
      }

      boolean isChiral = true;
      if( nInverted > 0 && ss.SingleMatch(tmpMol2)) isChiral = false;
      ss.delete();

      if( isChiral ) return;

      for(OEAtomBase orgAt : ringChiral)
      {  if( ! orgAt.GetBoolData(OEStruchk.ISStrainedBridgeHead)) continue;

         orgAt.SetBoolData(OEStruchk.STEREOClearTag, true);

         OEAtomBase tmpAt = getAtomByNum(copyAllStereoSpecified, orgAt);
         Atom.removeChiralInfo(tmpAt);
      }
   }


   /** get an atom from the molecule objects by the value stored in ATOMNumTag */
   static OEAtomBase getAtomByNum(OEMolBase mol, OEAtomBase srcMolAt)
   {  AtomIntPropertyFunctor atFunct = new AtomIntPropertyFunctor(
               OEStruchk.ATOMNumTag, srcMolAt.GetIntData(OEStruchk.ATOMNumTag));

      OEAtomBase tmpAt = mol.GetAtom(atFunct);
      atFunct.delete();

      return tmpAt;
   }


   private static OESubSearch buildSubSearch(OEGraphMol in)
   {  OEQMol qMol = new OEQMol(in);
      oechem.OEMDLClearParity(qMol);
      qMol.BuildExpressions(OEExprOpts.ExactAtoms|OEExprOpts.SymmetryClass,
                            OEExprOpts.ExactBonds|OEExprOpts.SymmetryClass);
      OESubSearch subSearch = new OESubSearch(qMol);
      qMol.delete();
      return subSearch;
   }
}


/**
 * Used to generate correct stereochemistry on strained  ring
 * systems such as adamantyl.
 *
 * @author albertgo
 *
 */
class StrainedCycleStereoGenerator
{  private final String smi;
   private final OESubSearch cycleSubSearch;
   private final int maxAtomIdx;
   private final OEAtomBaseVector qAtVec = new OEAtomBaseVector();
   private final OEAtomBaseVector tAtVec = new OEAtomBaseVector();

   StrainedCycleStereoGenerator(String smiles)
   {  this.smi = smiles;

      OEGraphMol cycMol = new OEGraphMol();
      oechem.OEParseSmiles(cycMol, smiles);
      oechem.OESuppressHydrogens(cycMol);

      this.maxAtomIdx = cycMol.GetMaxAtomIdx();
      OEQMol qMol = new OEQMol(cycMol);
      qMol.BuildExpressions(OEExprOpts.RingMember,
                            OEExprOpts.RingMember);

      this.cycleSubSearch = new OESubSearch(qMol);
      qMol.delete();
      cycMol.delete();
   }

   String getSmiles() { return smi;}

   /**
    * Assign stereochemistry to all bridgehead atoms in strained cycle.
    *
    * Also flag matched bridgehead atoms as OEStruchk.ISStrainedBridgeHead.
    */
   void assignStereo( OEMolBase mol, OEMolBase parentMol )
   {  OEAtomBase[] queryIdxToTargetMap = new OEAtomBase[maxAtomIdx];

      OEMatchBaseIter matIt = cycleSubSearch.Match(mol, true);
      while( matIt.hasNext() )   // for each ring system in mol
      {  OEMatchBase mat = matIt.next();

         // construct mapping from query to mol atoms
         OEMatchPairAtomIter atPIt = mat.GetAtoms();
         while( atPIt.hasNext() )
         {  OEMatchPairAtom atp = atPIt.next();
            OEAtomBase qAt = atp.getPattern();
            OEAtomBase tAt = atp.getTarget();
            atp.delete();

            queryIdxToTargetMap[qAt.GetIdx()] = tAt;
         }

         // if query had stereo set stereo on mol atom
         atPIt.ToFirst();
         while( atPIt.hasNext() )
         {  OEMatchPairAtom atp = atPIt.next();
            OEAtomBase qAt = atp.getPattern();
            if( ! qAt.HasStereoSpecified() )
            {  atp.delete();
               continue;
            }

            qAtVec.clear();
            tAtVec.clear();
            OEAtomBase tAt = atp.getTarget();
            atp.delete();

            OEAtomBaseIter qAtIt = qAt.GetAtoms();
            while( qAtIt.hasNext() )
            {  OEAtomBase qNeigh = qAtIt.next();
               qAtVec.add(qNeigh);
               tAtVec.add(queryIdxToTargetMap[qNeigh.GetIdx()]);
            }
            qAtIt.delete();

            // set stereo
            int qStereo = qAt.GetStereo(qAtVec, OEAtomStereo.Tetrahedral);
            tAt.SetStereo(tAtVec, OEAtomStereo.Tetrahedral, qStereo);
            tAt.DeleteData(OEStruchk.STEREOClearTag); // N in C[C@@]([O])([H])[C@]1([H])N2C[C@@](C=C)([H])[C@@H](CC2)C1

            OEAtomBase parentAt = FlagNonChiralStereoCenters.getAtomByNum(parentMol, tAt);
            parentAt.SetBoolData(OEStruchk.ISStrainedBridgeHead, true);
            tAt.SetBoolData(OEStruchk.ISStrainedBridgeHead, true);

            // N in C[C@@]([O])([H])[C@]1([H])N2C[C@@](C=C)([H])[C@@H](CC2)C1
            // was defined as not having stereo but now we know it does because it is
            // in a strained bridgehead position
            if( parentAt.DeleteData(OEStruchk.STEREOClearTag) )
            {  tAt.DeleteData(OEStruchk.STEREOClearTag);
               tAt.SetChiral(true);
               parentAt.SetChiral(true);
            }

         }

         atPIt.delete();
         mat.delete();
      }
      matIt.delete();
   }

   void delete()
   {  qAtVec.delete();
      tAtVec.delete();
      cycleSubSearch.delete();
   }
}
