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

import org.jdom.Element;

import com.aestel.utility.Message;
import com.aestel.utility.MessageList;
import com.genentech.oechem.tools.Atom;
import com.genentech.struchk.oeStruchk.OEStruchk.StructureFlag;

/**
 * Assume all stereo specification in molecule are absolute and assign stereo flag
 * based on the number of unspecified stereocenters.
 *
 * @author albertgo 2008
 *
 */
public class AssignStructureFlag extends AbstractStructureFlagCheck {

   private boolean structFlagIsCertain;


   public AssignStructureFlag(Element elem, FlagNonChiralStereoCenters stereoAtomFlagger) {
      super(elem, stereoAtomFlagger);
   }


   @Override
   public HydrogenMode getRequiredHydrogenMode() {
      return HydrogenMode.SUPRRESSED;
   }

   @Override
   public boolean checkStructure(OEGraphMol in, StructureFlag inStereo,
         MessageList msgs) {

      // make sure all special non-chiral atoms are flagged
      stereoAtomFlagger.checkStructure(in, null, msgs);

      boolean success = true;
      flag = StructureFlag.NOStereo;
      structFlagIsCertain = false;
      nChiralSpecified = 0;
      nChiralTotal = 0;
      nNonChiralStereoTotal = 0;
      nNonChiralStereoSpecified = 0;
      nStereoDBondSpecified = 0;
      nStereoDBondTotal = 0;
      int nAtoms = 0;

      // count atoms with stereo info and atoms which can have stereo info
      OEAtomBaseIter aIt = in.GetAtoms();
      while( aIt.hasNext() )  {
         OEAtomBase at = aIt.next();
         if(at.GetAtomicNum() != 0) // ignore "*" atoms
            nAtoms++;

//         // clear Nitrogen chirality
//         if( at.IsNitrogen() && at.IsChiral() && ! at.GetBoolData(OEStruchk.ISChiralNotRecognized) ) {
//            at.SetBoolData(OEStruchk.STEREOClearTag, true);
//
//            continue;
//         }
//
         // clear stereo on S+
         if( at.IsSulfur() && at.IsChiral() &&  at.GetFormalCharge() != 0 ) {
            at.SetBoolData(OEStruchk.STEREOClearTag, true);

            continue;
         }

         // Flagged as cleared by struCheck because of unstable atom
         if(at.GetBoolData(OEStruchk.STEREOClearTag)) {
            if( at.HasStereoSpecified() ) {
               Atom.removeChiralInfo(at);
               msgs.addMessage(new Message(String.format(
                     "Atom %s has no stereo and may not have wedge bonds.",
                     Atom.getAtomName(at)),
                     Message.Level.ERROR, null));
                  success = false;
            }
            continue;
         }


         if(at.IsChiral() || at.GetBoolData(OEStruchk.ISChiralNotRecognized)) {
            nChiralTotal++;
            if(at.HasStereoSpecified()) nChiralSpecified++;
         }

         if(at.GetBoolData(OEStruchk.NONChiralStereoAtomTag)) {
            nNonChiralStereoTotal++;
            if(at.HasStereoSpecified()) nNonChiralStereoSpecified++;
         }
      }
      aIt.delete();

      // count double bonds with stereo info and double bonds which can have stereo info
      OEBondBaseIter bIt = in.GetBonds();
      while(bIt.hasNext()) {
         OEBondBase bd = bIt.next();

         // Flagged as cleared by struCheck because unstable double bond
         if(bd.GetBoolData(OEStruchk.STEREOClearTag)) continue;

         if(bd.IsInRing() && oechem.OEBondGetSmallestRingSize(bd) <= 7) continue;

         if(bd.IsChiral()) {
            nStereoDBondTotal++;
            if(bd.HasStereoSpecified()) nStereoDBondSpecified++;
         }
      }
      bIt.delete();

      if(  nNonChiralStereoTotal != nNonChiralStereoSpecified
         ||nStereoDBondTotal     != nStereoDBondSpecified ) {
         // unspecified double bonds or sp3 centers
         // could also be SINGLEUnknownStereoIsomer, MIXTUREOfENantiomers
         flag = StructureFlag.MIXTUREOfDiastereomers;
      }

      // all (if any) double bonds and non chiral sp3 centers are specified

      else if(nChiralTotal == 0) {
         if(nStereoDBondTotal == 0 && nNonChiralStereoTotal == 0) {
            // must be
            flag = StructureFlag.NOStereo;
            structFlagIsCertain = true;

            if( nAtoms == 0 ) // special case NullStruct is usually an error checked elsewhere
               flag = StructureFlag.UNCERTAINStructure;

         } else {
            // no chiral atoms and other centers are specified, must be
            flag = StructureFlag.SINGLEStereoisomer;
            structFlagIsCertain = true;
         }
      }

      else if(nChiralTotal == 1) {
         if(nChiralTotal == nChiralSpecified) {
            // could also be SINGLEUnknownStereoIsomer or MIXTUREOfEnantiomers
            flag = StructureFlag.SINGLEStereoisomer;

         } else {
            // only one chiral center =>
            // could also be SINGLEUnknownStereoisomer.
            flag = StructureFlag.MIXTUREOfEnantiomers;
         }
      }

      else { //if(nTotalChiral > 1) {
         if(nChiralTotal == nChiralSpecified) {
            // could also be SINGLEUnknownStereoIsomer, MIXTUREOfEnantiomers
            // or MIXTUREOfDiastereomers
            flag = StructureFlag.SINGLEStereoisomer;

         } else {
            // only one chiral center =>
            // could also be a SINGLEUnknownStereoisomer or a MIXTUREOfEnantiomers.
            flag = StructureFlag.MIXTUREOfDiastereomers;
         }
      }

      // remove stereo centers for molecules which are not SINGLEStereoisomer
      if(flag != StructureFlag.SINGLEStereoisomer && flag != StructureFlag.NOStereo)
      {  removeChiralInfo(in);
         if(flag != StructureFlag.MIXTUREOfEnantiomers) {
            removeRingSP3Stereo(in);
            removeDBStereo(in);
         }
         oechem.OEMDLPerceiveBondStereo(in);
      }

      return success;
   }


   /* (non-Javadoc)
    * @see StructFlagAnalysisInterface#getStereoInfo()
    */
   @Override
   public StructureFlag getStructureFlag() {
      return flag;
   }


   @Override
   public String getDescription() {
      return "Assume all centers with stereo information are absolute and assign" +
            " the corresponding stereo information.";
   }


   public boolean structFlagIsCertain() {
      return structFlagIsCertain;
   }


   @Override
   public boolean hasStructureFlagError()
   {  return false;
   }
}
