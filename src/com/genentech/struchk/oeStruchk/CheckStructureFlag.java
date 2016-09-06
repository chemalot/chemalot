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

import java.util.ArrayList;
import java.util.List;

import openeye.oechem.*;

import org.jdom.Element;

import com.aestel.utility.Message;
import com.aestel.utility.MessageList;
import com.genentech.oechem.tools.Atom;
import com.genentech.struchk.oeStruchk.OEStruchk.StructureFlag;

/**
 * Check to see if all stereocenters have the correct stereo chemistry.
 *
 * This needs SuppressedHydrogens!
 *
 * The rules for checking the stereo flag depend on the number of specified and toatl
 * chiral atoms, stereo double bonds and cyclic sp3 centers:<br/>
 * <table border='1' style='border-collapse:collapse;text-align:center;font-size: 10pt;'>
 *  <tr>
 *   <td colspan=11 style='text-align:left;'>
 *      <b>Mapping of specified and unspecified stereo centers to valid gneStructureFlag</b><br/>
 *      N chiral: number of chiral atoms<br/>
 *      N db: number of asymmetric double bonds<br/>
 *      N non chiral stereo number of sp3 atoms which are not chiral but bear stereo info eg. dimethyl-cyclo-hexan</td>
 *  </tr>
 *  <tr>
 *   <td colspan=2>N chiral</td>
 *   <td colspan=2>N db</td>
 *   <td colspan=2>N non chiral stereo</td>
 *   <td rowspan=2>No Stereo</td>
 *
 *   <td rowspan=2>Single Stereoisomer</td>
 *   <td rowspan=2>Single Unknown SI</td>
 *   <td rowspan=2>Mix Enantiomer</td>
 *   <td rowspan=2>Mix Diastereomers</td>
 *  </tr>
 *  <tr>
 *   <td>spc</td><td>Nsp</td><td>spc</td><td>nsp</td><td>spc</td><td>nsp</td>
 *  </tr>
 *
 *  <tr>
 *   <td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td>
 *   <td>OK</td><td>Error B</td><td>Error B</td><td>Error B</td><td>Error B</td>
 *  </tr>
 *  <tr>
 *   <td>1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td>
 *   <td>Error A</td><td>OK</td><td>OK</td><td>OK</td><td>Error G</td>
 *  </tr>
 *  <tr>
 *   <td>&gt;=2</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td>
 *   <td>Error A</td><td>OK</td><td>OK</td><td>OK</td><td>OK</td>
 *  </tr>
 *  <tr>
 *   <td>0</td><td>1</td><td>0</td><td>0</td><td>0</td><td>0</td>
 *   <td>Error A</td><td>Error C</td><td>OK</td><td>OK</td><td>Error G</td>
 *  </tr>
 *  <tr>
 *   <td>0</td><td>&gt;=2</td><td>0</td><td>0</td><td>0</td><td>0</td>
 *   <td>Error A</td><td>Error C</td><td>OK</td><td>OK</td><td>OK</td>
 *  </tr>
 *  <tr>
 *   <td>0</td><td>0</td><td>&gt;=1</td><td>0</td><td>0</td><td>0</td>
 *   <td>Error A</td><td>OK</td><td>Error D</td><td>Error D+E</td><td>Error D</td>
 *  </tr>
 *  <tr>
 *   <td>0</td><td>0</td><td>0</td><td>&gt;=1</td><td>0</td><td>0</td>
 *   <td>Error A</td><td>Error C</td><td>OK</td><td>Error E+F</td><td>OK</td>
 *  </tr>
 *  <tr>
 *   <td>&gt;=1</td><td>&gt;=1</td><td>0</td><td>0</td><td>0</td><td>0</td>
 *   <td>Error A</td><td>Error C</td><td>OK</td><td>OK</td><td>OK</td>
 *  </tr>
 *  <tr>
 *   <td>1</td><td>0</td><td>&gt;=1</td><td>0</td><td>0</td><td>0</td>
 *   <td>Error A</td><td>OK</td><td>OK</td><td>OK</td><td>OK</td>
 *  </tr>
 *  <tr>
 *   <td>&gt;=2</td><td>0</td><td>&gt;=1</td><td>0</td><td>0</td><td>0</td>
 *   <td>Error A</td><td>OK</td><td>OK</td><td>OK</td><td>OK</td>
 *  </tr>
 *  <tr>
 *   <td>&gt;=1</td><td>0</td><td>0</td><td>&gt;=1</td><td>0</td><td>0</td>
 *   <td>Error A</td><td>Error C</td><td>OK</td><td>Error F</td><td>OK</td>
 *  </tr>
 *  <tr>
 *   <td>0</td><td>1</td><td>&gt;=1</td><td>0</td><td>0</td><td>0</td>
 *   <td>Error A</td><td>Error C</td><td>OK</td><td>OK</td><td>Error G</td>
 *  </tr>
 *  <tr>
 *   <td>0</td><td>&gt;2</td><td>&gt;=1</td><td>0</td><td>0</td><td>0</td>
 *   <td>Error A</td><td>Error C</td><td>OK</td><td>OK</td><td>OK</td>
 *  </tr>
 *  <tr>
 *   <td>0</td><td>&gt;=1</td><td>0</td><td>&gt;=1</td><td>0</td><td>0</td>
 *   <td>Error A</td><td>Error C</td><td>OK</td><td>Error F</td><td>OK</td>
 *  </tr>
 *  <tr>
 *   <td>0</td><td>0</td><td>&gt;=1</td><td>&gt;=1</td><td>0</td><td>0</td>
 *   <td>Error A</td><td>Error C</td><td>OK</td><td>Error E</td><td>OK</td>
 *  </tr>
 *  <tr>
 *   <td colspan=12>&nbsp;</td>
 *  </tr>
 *  <tr>
 *   <td>0</td><td>0</td><td>0</td><td>0</td><td>&gt;=2</td><td>0</td>
 *   <td>Error A</td><td>OK</td><td>Error D</td><td>Error D+E</td><td>Error D</td>
 *  </tr>
 *  <tr>
 *   <td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>&gt;=2</td>
 *   <td>Error A</td><td>Error C</td><td>OK</td><td>Error E+F</td><td>OK</td>
 *  </tr>
 *  <tr>
 *   <td>1</td><td>0</td><td>0</td><td>0</td><td>&gt;=2</td><td>0</td>
 *   <td>Error A</td><td>OK</td><td>OK</td><td>OK</td><td>Error G</td>
 *  </tr>
 *  <tr>
 *   <td>&gt;=2</td><td>0</td><td>0</td><td>0</td><td>&gt;=2</td><td>0</td>
 *   <td>Error A</td><td>OK</td><td>OK</td><td>OK</td><td>OK</td>
 *  </tr>
 *  <tr>
 *   <td>&gt;=1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>&gt;=2</td>
 *   <td>Error A</td><td>Error C</td><td>OK</td><td>Error F</td><td>OK</td>
 *  </tr>
 *  <tr>
 *   <td>0</td><td>1</td><td>0</td><td>0</td><td>&gt;=2</td><td>0</td>
 *   <td>Error A</td><td>Error C</td><td>OK</td><td>OK</td><td>Error G</td>
 *  </tr>
 *  <tr>
 *   <td>0</td><td>&gt;2</td><td>0</td><td>0</td><td>&gt;=2</td><td>0</td>
 *   <td>Error A</td><td>Error C</td><td>OK</td><td>OK</td><td>OK</td>
 *  </tr>
 *  <tr>
 *   <td>0</td><td>&gt;=1</td><td>0</td><td>0</td><td>0</td><td>&gt;=2</td>
 *   <td>Error A</td><td>Error C</td><td>OK</td><td>Error F</td><td>OK</td>
 *  </tr>
 *  <tr>
 *   <td>0</td><td>0</td><td>0</td><td>0</td><td>&gt;=2</td><td>&gt;=2</td>
 *   <td>Error A</td><td>Error C</td><td>OK</td><td>Error E</td><td>OK</td>
 *  </tr>
 *  <tr>
 *   <td colspan=11 style='text-align:left;'>
 *    Error A: Mol has stereo centers may not be "No Stereo"<br>
 *    Error B: Mol has no stereo centers and must be "No Stereo"<br/>
 *    Error C: Mol has unspecified chiral centers, double bonds  or cyclic sp3 atoms.<br/>
 *    Error D: Mol has full stereo information and must be "Single Stereoisomer"<br/>
 *    Error E: Mol has no chiral atoms can not be "Mixture of Enantiomers"<br/>
 *    Error F: Mol has unspecified double bonds and/or cyclic sp3  atoms and can not be a "Mixture of Enantiomers".<br/>
 *    Error G: Mol has only one chiral center and is therefore a "Mixture of Enantiomers"</td>
 *  </tr>
 * </table>
 *
 * @author albertgo 2008
 *
 */
public class CheckStructureFlag extends AbstractStructureFlagCheck {

   private boolean hasFlagError = false;

   public CheckStructureFlag(Element elem, FlagNonChiralStereoCenters stereoAtomFlagger) {
      super(elem, stereoAtomFlagger);
   }


   @Override
   public HydrogenMode getRequiredHydrogenMode() {
      return HydrogenMode.SUPRRESSED;
   }


   @Override
   @SuppressWarnings("unchecked")
   public boolean checkStructure(OEGraphMol in, StructureFlag inFlag,
                                 MessageList msgs) {
      // we are currently not using the OEMDLHasParity function (the MDL chiral flag)

      assert inFlag != null;

      nChiralSpecified = 0;
      nChiralTotal = 0;
      nNonChiralStereoTotal = 0;
      nNonChiralStereoSpecified = 0;
      nStereoDBondSpecified = 0;
      nStereoDBondTotal = 0;

      // make sure all special non-chiral atoms are flagged
      stereoAtomFlagger.checkStructure(in, inFlag, msgs);

      flag = inFlag;
      boolean success = true;
      hasFlagError = false;

      // count atoms with stereo info and atoms which can have stereo info
      int nAtoms = 0;
      OEAtomBaseIter aIt = in.GetAtoms();
      List<OEAtomBase> chiralAtoms = new ArrayList<OEAtomBase>(in.GetMaxAtomIdx());
      List<OEAtomBase> nonChiralAtoms = new ArrayList<OEAtomBase>(in.GetMaxAtomIdx());

      while( aIt.hasNext() )  {
         OEAtomBase at = aIt.next();
         if(at.GetAtomicNum() != 0) // ignore "*" atoms
            nAtoms++;

//         // clear Nitrogen chirality
//         if( at.IsNitrogen() && at.IsChiral() && ! at.GetBoolData(OEStruchk.ISChiralNotRecognized) ) {
//            at.SetBoolData(OEStruchk.STEREOClearTag, true);
//            continue;
//         }

         // clear stereo on S+
         if( at.IsSulfur() && at.IsChiral() &&  at.GetFormalCharge() != 0 ) {
            at.SetBoolData(OEStruchk.STEREOClearTag, true);

            continue;
         }

         // Flagged as cleared by struCheck because of unstable atom
         if(at.GetBoolData(OEStruchk.STEREOClearTag)) {
            if( at.HasStereoSpecified() && ! at.GetBoolData(OEStruchk.ATROPIsomericCenter)) {
               msgs.addMessage(new Message(String.format(
                     "Atom %s has no stereo and may not have wedge bonds.",
                     Atom.getAtomName(at)),
                     Message.Level.ERROR, null));
                  success = false;
               Atom.removeChiralInfo(at);
            }
            continue;
         }

         if(at.IsChiral() || at.GetBoolData(OEStruchk.ISChiralNotRecognized)) {
            nChiralTotal++;
            chiralAtoms.add(at);
            if(at.HasStereoSpecified()) nChiralSpecified++;
         }

         if(at.GetBoolData(OEStruchk.NONChiralStereoAtomTag)) {
            nNonChiralStereoTotal++;
            nonChiralAtoms.add(at);
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

      /* logic reproducing rules described in class javadoc */
      if(nChiralTotal == 0 && nStereoDBondTotal == 0 && nNonChiralStereoTotal == 0) {
         // special case: NullStruct with nAtoms == 0  is checked else where
         if(nAtoms > 0 && inFlag != StructureFlag.NOStereo && inFlag != StructureFlag.UNCERTAINStructure) {
            msgs.addMessage(new Message("Molecule has no stereo centers and must be flagged as such.",
                  Message.Level.ERROR, null));
            success = false;
            hasFlagError = true;
         } else if( nAtoms == 0 && inFlag != StructureFlag.UNCERTAINStructure ) {
            msgs.addMessage(new Message("Null structure can only be an 'Uncertain Strucutre'.",
                  Message.Level.ERROR, null));
            success = false;
         }
      } else if(inFlag == StructureFlag.NOStereo) {
         if(nChiralTotal != 0 || nNonChiralStereoTotal != 0 || nStereoDBondTotal != 0) {
               msgs.addMessage(new Message(String.format(
                  "Molecule has %d stereo centers %sand %d stereo double bonds and can not be specified as 'No Stereo'!",
                  nChiralTotal + nNonChiralStereoTotal, getAtomList(chiralAtoms, nonChiralAtoms),
                  nStereoDBondTotal), Message.Level.ERROR, null));
               success = false;
               hasFlagError = true;
         }
      } else if(inFlag == StructureFlag.SINGLEStereoisomer) {
         if(nChiralTotal != nChiralSpecified
            || nNonChiralStereoTotal != nNonChiralStereoSpecified
            || nStereoDBondTotal != nStereoDBondSpecified) {
            msgs.addMessage(new Message(String.format(
                  "This Single Stereoisomer has %d stereo centers %sand %d stereo double bonds which all must be specified!",
                  nChiralTotal + nNonChiralStereoTotal, getAtomList(chiralAtoms, nonChiralAtoms),
                  nStereoDBondTotal), Message.Level.ERROR, null));
               success = false;
               hasFlagError = true;
         }
      } else if(inFlag == StructureFlag.SINGLEUnknownStereoIsomer) {
         if(nChiralTotal == 0
               && nNonChiralStereoTotal == nNonChiralStereoSpecified
               && nStereoDBondTotal == nStereoDBondSpecified) {
// users may draw SUS with or without stereo specification on any atoms and bonds
// eg. a user might have separated dimethylcyclohexyl but not know if it is cis or trans
//            msgs.addMessage(new Message(String.format(
//                  "This Molecule has all sp3 and double bonds specified and should be a 'Single Stereoisomer'!",
//                  nNonChiralStereo, nBondChiral), Message.Level.ERROR, null));
//            success = false;
//            hasFlagError = true;
         }
      } else if(inFlag == StructureFlag.MIXTUREOfEnantiomers) {
         if( nChiralTotal == 0 ) {
            msgs.addMessage(new Message(
                  "Molecule has no chiral centers and can not be specified as 'Mixture of Enantiomer'!",
                  Message.Level.ERROR, null));
            success = false;
            hasFlagError = true;
         }
         else if(   nNonChiralStereoTotal != nNonChiralStereoSpecified
            || nStereoDBondTotal != nStereoDBondSpecified) {
            msgs.addMessage(new Message(String.format(
               "Molecule has %d sp3 ring atoms %sand %d double bonds which all must have stereo chemistry!",
               nNonChiralStereoTotal, getAtomList(nonChiralAtoms),
               nStereoDBondTotal), Message.Level.ERROR, null));
            success = false;
            hasFlagError = true;
         }
      } else if(inFlag == StructureFlag.MIXTUREOfDiastereomers) {
         if(nChiralTotal == 0
               && nNonChiralStereoTotal == nNonChiralStereoSpecified
               && nStereoDBondTotal == nStereoDBondSpecified
               && nNonChiralStereoSpecified/2 + nStereoDBondTotal == 1) {
            // an MD with just one center of stereochemistry must be drawn flat
            // as that center has to be a mixture
            // If there are more than one center then one can be known absolutely
            // while the other is a mixture or multiple can be mixtures
            msgs.addMessage(new Message(String.format(
                  "This Molecule has %d stereodouble bonds and %d cis/trans rings. This has to be drawn without stereo!",
                  nStereoDBondTotal, nNonChiralStereoTotal/2), Message.Level.ERROR, null));
            success = false;
            hasFlagError = true;
         } else if(nChiralTotal == 1
               && nNonChiralStereoTotal == nNonChiralStereoSpecified
               && nStereoDBondTotal == nStereoDBondSpecified) {
            msgs.addMessage(new Message(String.format(
                  "This Molecule has one chiral atom and all sp3 and double bonds specified and should be a 'Mixture Of Enantiomers'!",
                  nNonChiralStereoTotal, nStereoDBondTotal), Message.Level.ERROR, null));
            success = false;
            hasFlagError = true;
         }
      }

      // If compound not pure stereoisomer remove stereo info from bonds and atoms
      if(inFlag != StructureFlag.SINGLEStereoisomer)
      {  removeChiralInfo(in);
         if(inFlag != StructureFlag.MIXTUREOfEnantiomers) {
            removeRingSP3Stereo(in);
            removeDBStereo(in);
         }
         oechem.OEMDLPerceiveBondStereo(in);
      }

      return success;
   }

   private String getAtomList(List<OEAtomBase>... atomLists) {
      StringBuilder sb = new StringBuilder(50);
      for(List<OEAtomBase>alist : atomLists) {
         for(OEAtomBase at : alist) {
            sb.append(Atom.getAtomName(at)).append(", ");
         }
      }

      if( sb.length() > 0 ) sb.setLength(sb.length()-2);
      if( sb.length() > 0 ) sb.insert(0, "(").append(") ");
      return sb.toString();
   }


   @Override
   public boolean hasStructureFlagError()
   {  return hasFlagError;
   }


   @Override
   public String getDescription() {
      return "Check that StructureFlag is compatible with the stereo specifications" +
             " in the molecule.";
   }
}
