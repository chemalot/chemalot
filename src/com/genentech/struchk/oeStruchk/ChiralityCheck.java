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

import java.util.HashSet;

import org.jdom.Element;

import com.aestel.utility.Message;
import com.aestel.utility.MessageList;
import com.genentech.oechem.tools.Atom;
import com.genentech.oechem.tools.BondHasDataFunctor;
import com.genentech.oechem.tools.BondOrderFunctor;
import com.genentech.struchk.oeStruchk.OEStruchk.StructureFlag;

import openeye.oechem.OEAndBond;
import openeye.oechem.OEAtomBase;
import openeye.oechem.OEAtomBaseIter;
import openeye.oechem.OEBondBase;
import openeye.oechem.OEBondBaseIter;
import openeye.oechem.OEGraphMol;
import openeye.oechem.OEProperty;
import openeye.oechem.oechem;

/**
 * Structure checking rules which flags an input compound has wedge bonds which
 * define an invalid or not existing stereo chemistry.
 *
 * @author albertgo
 *
 */
public class ChiralityCheck extends AbstractStructureCheck{
   private final BondHasDataFunctor bdDataFctr;
   private final OEAndBond checkBondStereoFctr;
   private final BondOrderFunctor bdOrderFctr;
   private final HashSet<Integer> badAtoms = new HashSet<Integer>(20);

   public ChiralityCheck(Element tElem) {
      super(tElem);

      if(! checkExample())
         throw new Error( String.format("Example %s did not have invalid atom stereo.",
                                        getExampleInput()));

      bdDataFctr = new BondHasDataFunctor(OEProperty.BondStereo);
      bdOrderFctr = new BondOrderFunctor(1);
      checkBondStereoFctr = new OEAndBond(bdDataFctr, bdOrderFctr);
   }


   @Override
   public boolean checkStructure(OEGraphMol in, StructureFlag inStereo,
                                 MessageList msgs) {
      boolean hasError = false;
      badAtoms.clear();

      // check for bonds which had wedges but are not chiral and are not exocyclic singles
      hasError = checkForAtropIsomersOrMetals(in, msgs, hasError);

      checkForUnsupportedChirality(in, msgs);
      checkForNonChiralAtoms(in, msgs);

      if(! hasError && oechem.OEMDLHasIncorrectBondStereo(in)) {
         msgs.addMessage(new Message(
               "Stereochemistry is specified inconsistently.", Message.Level.ERROR, null));
         hasError = true;
      }

      if(hasError) return false;

      return true;
   }


   // check for bonds which had wedges but are not chiral and are not exocyclic singles
   private boolean checkForAtropIsomersOrMetals(OEGraphMol in, MessageList msgs,
            boolean hasError) {

      oechem.OEFindRingAtomsAndBonds(in);
      OEBondBaseIter bdIt = in.GetBonds(checkBondStereoFctr);
      while(bdIt.hasNext()) {
         OEBondBase bd = bdIt.next();
         OEAtomBase at = bd.GetBgn();
         
         if( at.IsChiral() ) continue;
         if( at.IsMetal() && at.GetDegree() > 4)
         {  at.SetBoolData(OEStruchk.CHIRALMetal, true);
            continue;
         }
         
         // we rely on openeye to remove stereo chemistry from atoms which do not have any
         if( at.GetBoolData(OEStruchk.NONChiralStereoAtomTag) ) continue; // exocyclic single bonds are included

         boolean isAtropIsomer = isAtropIsomericWedge(at, bd);
         if( isAtropIsomer ) {
            if( at.HasStereoSpecified() ) Atom.removeChiralInfo(at);
            at.SetChiral(false);

            at.SetBoolData(OEStruchk.STEREOClearTag, true);
            continue;
         }
         if( at.GetBoolData(OEStruchk.STEREOClearTag) ) continue;


         int atIdx = at.GetIdx();
         if( badAtoms.contains(atIdx) ) continue;

         hasError = true;
         badAtoms.add(atIdx);
         if( at.HasStereoSpecified() ) Atom.removeChiralInfo(at);

         msgs.addMessage(new Message(String.format("Atom %s is not chiral but has wedge bonds.",
               Atom.getAtomName(at)),
               Message.Level.ERROR,at));

      }
      bdIt.delete();
      return hasError;
   }

   /**
    * Check to see if the wedge bd is a valid assuming that atropAt is atropisomeric.
    *
    *  - atropAt must have 3 bonds (explicit to enfoce size restriction)
    *  - most have a single bond that connects it to annother atom with 3 bonds
    *  - That connecting bond may not be the wedgeBd and is the atropisomeric center
    *
    * @param atropAt atom with narrow side of wedge this would be one side of the atropisomeric center
    * @param wedgeBd on which wedge bond was drawn
    *
    * @return true if this atropAt could be part of an atropisomeric center
    *
    *
    */
   private boolean isAtropIsomericWedge(OEAtomBase atropAt, OEBondBase wedgeBd) {
      if( atropAt.GetImplicitHCount() == 0 && atropAt.GetExplicitDegree() != 3 ) return false;

      boolean isAtropIsomeric = false;

      // Atropisomer must be between two sp2 atoms with 3 substituents
      if( atropAt.GetHvyDegree() != 3 || atropAt.GetImplicitHCount() != 0)
         return isAtropIsomeric;

      // check for AtropIsomer
      OEBondBaseIter bit2 = atropAt.GetBonds();
      while(bit2.hasNext()) {
         OEBondBase b2 = bit2.next();
         if( b2.IsAromatic() || b2.GetIntType() != 1 )
            continue;  // we are looking for the single aliphatic bond connecting the two atoms

         int b2Index = b2.GetIdx();
         if( b2Index == wedgeBd.GetIdx()) continue;

         int smalRingSize = oechem.OEBondGetSmallestRingSize(b2);
         if( smalRingSize > 0 && smalRingSize < 7 ) {
            // atropisomeric bond may not be in a small ring (<7) because this would
            // cause the atoms to be close to planar thus reducing the transition state energy
            continue;
         }

         OEAtomBase at2 = b2.GetNbr(atropAt);
         // Atropisomer must be between two sp2 atoms with 3 substituents
         if( at2.GetImplicitHCount() == 0 && at2.GetHvyDegree() == 3 ) {

            isAtropIsomeric = true;

            // look for a wedge bond on the other side of the atropIsomeirc bond (b2)
            // only one wedge is allowed specifying the stereo of an atropisomer
            OEBondBaseIter bit3 = at2.GetBonds();
            while(bit3.hasNext()) {
               OEBondBase b3 = bit3.next();
               if( b3.GetIdx() == b2Index ) continue;

               if( b3.GetIntType() == 1 && b3.HasData(OEProperty.BondStereo) ) {
                  OEAtomBase at3 = b3.GetBgn();
                  if( at3.GetIdx() == at3.GetIdx() ) {
                     // other side of atrop bond has also a wedge bond
                     // we consider this an error since it makes identifying
                     // the 3d geometry difficult
                     // TODO: Add Message saying the two wedge bonds are not allowed across atropiisomeric centers
                     //       This should not be a problem even when we have nighboring atrobisomeric bonds.
                     isAtropIsomeric = false;
                  }
               }
            }
            bit3.delete();

            if( isAtropIsomeric ) {
               // valid connectivity for atropIsomer
               atropAt.SetBoolData(OEStruchk.ATROPIsomericCenter, true);
               at2.SetBoolData(OEStruchk.ATROPIsomericCenter, true);
               b2.SetBoolData(OEStruchk.ATROPIsomericCenter, true);

               break;
            }
         }
      }
      bit2.delete();

      return isAtropIsomeric;
   }

   private void checkForUnsupportedChirality(OEGraphMol in, MessageList msgs) {

      OEAtomBaseIter aIt = in.GetAtoms();
      while( aIt.hasNext() ) {
         OEAtomBase at = aIt.next();

         if(at.HasStereoSpecified() && ! at.GetBoolData(OEStruchk.ATROPIsomericCenter)) {
            if( ! FlagNonChiralStereoCenters.canBeChiral(at) ) {
               Atom.removeChiralInfo(at);

               msgs.addMessage(new Message(String.format("Chirality on atom %s not supported.",
                  Atom.getAtomName(at)),
                  Message.Level.ERROR,at));
               badAtoms.add(at.GetIdx());

            } else if( at.IsAromatic() ) {
               Atom.removeChiralInfo(at);

               msgs.addMessage(new Message(String.format("Chirality on aromatic atom %s not supported.",
                  Atom.getAtomName(at)),
                  Message.Level.ERROR,at));
               badAtoms.add(at.GetIdx());
            }
         }
      }
      aIt.delete();
   }

   private boolean checkForNonChiralAtoms(OEGraphMol in, MessageList msgs) {
      boolean hasError = false;

      OEAtomBaseIter aIt = in.GetAtoms();
      while( aIt.hasNext() ) {
         OEAtomBase at = aIt.next();

         if(at.HasStereoSpecified() && ! at.IsChiral() ) {

            // checkAChiralStereoAtom necessary for C[C@@H]1CC[C@H](CC1)C
            if(   ! at.GetBoolData(OEStruchk.NONChiralStereoAtomTag)
               && ! at.GetBoolData(OEStruchk.ATROPIsomericCenter)
               && ! badAtoms.contains(at.GetIdx())) {
               hasError = true;

               Atom.removeChiralInfo(at);

               msgs.addMessage(new Message(String.format("Atom %s is not chiral but has wedge bonds.",
                     Atom.getAtomName(at)),
                     Message.Level.ERROR,at));

               badAtoms.add(at.GetIdx());
            }
         }
      }
      aIt.delete();
      return hasError;
   }


   @Override
   public void delete() {
      bdDataFctr.delete();
      bdOrderFctr.delete();
      checkBondStereoFctr.delete();
   }
}
