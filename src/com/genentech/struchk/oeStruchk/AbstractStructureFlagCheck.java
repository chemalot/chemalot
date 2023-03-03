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

import org.jdom.Element;

import openeye.oechem.*;

import com.aestel.utility.MessageList;
import com.genentech.oechem.tools.Atom;
import com.genentech.oechem.tools.Bond;
import com.genentech.struchk.oeStruchk.OEStruchk.StructureFlag;

abstract class AbstractStructureFlagCheck implements StructFlagAnalysisInterface
{  protected final String checkName;
   protected final FlagNonChiralStereoCenters stereoAtomFlagger;
   protected final OEStruchk parentChecker;
   protected StructureFlag flag;
   protected int nStereoDBondSpecified;
   protected int nStereoDBondTotal;
   protected int nNonChiralStereoSpecified;
   protected int nNonChiralStereoTotal;
   protected int nChiralSpecified;
   protected int nChiralTotal;
   protected int nChiralNonTetrahedral;
   protected int nChiralNonTetrahedralSpecified;
   private OEGraphMol keeperSubstance;
   private MessageList keeperMsgs;

   protected AbstractStructureFlagCheck(OEStruchk checker, Element elem, FlagNonChiralStereoCenters stereoAtomFlagger) {
      checkName = elem.getName();
      parentChecker = checker;

      if( stereoAtomFlagger == null )
         throw new Error("CheckStructureFlag needs flagNonChiralAtoms");

      this.stereoAtomFlagger = stereoAtomFlagger;
   }

   /* remove sterochemistry from chiral atoms, leave nonchiral sp3 atoms alone */
   protected void removeChiralInfo(OEGraphMol in)
   {  int ringSysByAtom[] = null;
      OEAtomBaseIter aIt = in.GetAtoms();
      while( aIt.hasNext() )  {
         OEAtomBase at = aIt.next();

         if( at.HasData( OEStruchk.ATOMExpHX )) {
            at.DeleteData(OEStruchk.ATOMExpHX );
            at.DeleteData(OEStruchk.ATOMExpHY );
            at.DeleteData(OEStruchk.ATOMExpHZ );
         }

         if(!at.HasStereoSpecified() ) continue;
         if( ! at.IsChiral() )
         {  // some Nonchiral sp3 center depend on chiral center to be valid
            // eg. in: [O][C@@H]1C[C@H](C)[C@H](C)C1
            // their stereo needs to be removed as well
            // the condition is that the ring system needs to have more than one
            // NONChiralStereoAtom

            if( at.IsInRing() && at.GetBoolData(OEStruchk.NONChiralStereoAtomTag ) )
            {  if( ringSysByAtom == null )
               {  ringSysByAtom = new int[in.GetMaxAtomIdx()];
                  oechem.OEDetermineRingSystems(in, ringSysByAtom);
               }

               int nNChiralInRing = 0;
               int at1RingSystem = ringSysByAtom[at.GetIdx()];
               OEAtomBaseIter atIt2 = in.GetAtoms();
               while( atIt2.hasNext() )
               {  OEAtomBase at2 = atIt2.next();
                  if( ringSysByAtom[at2.GetIdx()] == at1RingSystem
                      && at2.GetBoolData(OEStruchk.NONChiralStereoAtomTag))
                     nNChiralInRing++;
               }
               atIt2.delete();

               if( nNChiralInRing > 1)
                  continue;   // this is a valid NONChiral sp3 atom as in CC1CC(C)C1
            }
         }

         Atom.removeChiralInfo(at);
      }
      aIt.delete();

      //nChiralSpecified = 0;
   }

   protected void removeRingSP3Stereo(OEGraphMol in)
   {  OEAtomBaseIter aIt = in.GetAtoms();
      while( aIt.hasNext() )  {
         OEAtomBase at = aIt.next();

         if(   ! at.HasStereoSpecified()
            || ! at.GetBoolData(OEStruchk.NONChiralStereoAtomTag) ) continue;

         Atom.removeChiralInfo(at);
      }
      aIt.delete();
      //nNonChiralStereoSpecified = 0;
   }

   protected void removeDBStereo(OEGraphMol in) {
      OEBondBaseIter bIt = in.GetBonds();
      while( bIt.hasNext() ) {
         OEBondBase bd = bIt.next();
         if(bd.GetOrder() != 2 || ! bd.HasStereoSpecified()) continue;

         Bond.removeDBStereo(bd);
      }
      bIt.delete();
      //nStereoDBondSpecified = 0;
   }

   /** return true if at matches O=*-[O;H] */
   static boolean hasOOHNeighbor(OEAtomBase at) {
      boolean hasOH = false;
      boolean hasdoubleO = false;
      
      OEBondBaseIter bdIt = at.GetBonds();
      while(bdIt.hasNext()) {
         OEBondBase bd = bdIt.next();
         OEAtomBase oAt = bd.GetNbr(at);
         
         if( ! oAt.IsOxygen() ) continue;

         int bdType = bd.GetIntType();
         if(  bdType == 1) {
            if( oAt.GetTotalHCount() > 0 ) {
               hasOH = true;
               if( hasdoubleO ) break;
            }
         } else if( bdType == 2 ) {
            hasdoubleO = true;
            if( hasOH ) break;
         }
      }   
      bdIt.delete();
      return hasOH && hasdoubleO;
    }



   @Override
   public int getNChiral()
   {  return nChiralTotal;
   }

   @Override
   public int getNChiralSpecified()
   {  return nChiralSpecified;
   }

   @Override
   public int getNChiralNonTetrahedral()
   {  return nChiralNonTetrahedral;
   }

   @Override
   public int getNChiralNonTetrahedralSpecified()
   {  return nChiralNonTetrahedralSpecified;
   }

   /**
    * @return total number of non chiral sp3 centers which can have wedge bonds
    *         eg. 1,4 dimethyl cyclohexan.
    */
   @Override
   public int getNNonChiralStereo()
   {  return nNonChiralStereoTotal;
   }

   /**
    * @return number of non chiral sp3 centers which have wedge bonds
    *         eg. trans 1,4 dimethyl cyclo hexan.
    */
   @Override
   public int getNNonChiralStereoSpecified()
   {  return nNonChiralStereoSpecified;
   }

   @Override
   public int getNStereoDBond()
   {  return nStereoDBondTotal;
   }

   @Override
   public int getNStereoDBondSpecified()
   {  return nStereoDBondSpecified;
   }

   @Override
   public StructureFlag getStructureFlag()
   {  return flag;
   }

   @Override
   public String getCheckName() {
      return checkName;
   }


   @Override
   public void reset() {
      flag = StructureFlag.NOStereo;
      nChiralNonTetrahedral = 0;
      nChiralNonTetrahedralSpecified = 0;
      nChiralSpecified = 0;
      nChiralTotal = 0;
      nNonChiralStereoTotal = 0;
      nNonChiralStereoSpecified = 0;
      nStereoDBondSpecified = 0;
   }

   @Override
   public String getDescriptionHTML() {
      return getDescription();
   }


   @Override
   public boolean checkExample() {
      return true;
   }

   @Override
   public String getExampleInput() {
      return null;
   }

   /** StructureKeeperInterface **/

   void keepSubstance(OEGraphMol in, MessageList msgs) {
      if( keeperSubstance != null)
      {  keeperSubstance.delete();
         keeperSubstance = null;
      }
      keeperSubstance = new OEGraphMol(in); // keep a copy so nobody can change it

      oechem.OESuppressHydrogens(keeperSubstance,false,false,true);
      oechem.OEPerceiveChiral(keeperSubstance);

      this.keeperMsgs = new MessageList(msgs);
   }

   @Override
   public OEGraphMol getMolecule() {
      return keeperSubstance;
   }

   @Override
   public String getKeeperName() {
      return STEREONormalizedKeeper;
   }


   @Override
   public MessageList getStructureMessages() {
      return keeperMsgs;
   }

   @Override
   public void delete() {
      if( keeperSubstance != null) keeperSubstance.delete();
      keeperSubstance = null;
   }
}
