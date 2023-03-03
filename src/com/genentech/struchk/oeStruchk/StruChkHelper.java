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

import java.util.EnumSet;
import java.util.Set;

public class StruChkHelper {
   public static enum CHECKType {
      atomLabelCheck,
      thickBondCheck,
      flagNonChiralAtoms,
      wigglyBondCheck,
      closeAtomCheck,
      twoDCheck,
      cleanReactionCenter,
      chiralCheck,
      clearAmidine1Stereo,
      clearAmidine2Stereo,
      removeHydrogens,
      checkDoubleBond,
      checkAtomTypes,
      separateAlkaliMetals,
      chargeAlkaliMetals,
      separateEarthAlkaliMetals,
      chargeEarthAlkaliMetals,
      amoniumHalides,
      protonateOS,
      protonateHalogen,
      deprotonateN,
      nitro,
      nitricAcid,
      noxide,
      azide,
      thiocyanate,
      sulforousAcid,
      perchloricAcid,
      sulfoxide,
      ketoEnol,
      checkAromaticTautomer,
      keepSubstance,
      dbComponentNormalizer,
      testComponentNormalizer,
      valenceCheck,
      keepAllStereo,
      assignStructFlag,
      checkStructFlag,
      keepParent,
      tautomerize,
      keepTautomer
   }

   /**
    * Enum of predefined structure checker configurations.
    */
   public static enum CHECKConfig {
      /** A CHECKConfig which assumes that all stereo centers
       * have been assigned absolutely and infers the StructFlag flag from
       * the stereo centers.
       *
       * It uses {@link AssignStructureFlag} and {@link ComponentNormalizer}.
       */
      ASSIGNStructFlag(StruChkHelper.ASSIGNStructFlagSet),

      /** A CHECKConfig which expects a StructFlag flag and will do its best to
       * validate that the StructFlag flag with the stereochemistry specified.
       *
       * It uses {@link CheckStructureFlag} and {@link ComponentNormalizer}.
       */
      CHECKStructFlag(StruChkHelper.CHECKStructFlagSet),


      /** These are flags for test configurations which should not be used in
       * production.
       *
       * It uses {@link CheckStructureFlag} and {@link ComponentNormalizer} but
       * it reads the salt info from the xml file.
       */
      TESTCheckStructFlag(StruChkHelper.TESTCheckStructFlagSet),

      /** These are flags for test configurations which should not be used in
       * production.
       *
       * It uses {@link AssignStructureFlag} and {@link ComponentNormalizer} but
       * it reads the salt info from the xml file.
       */
      TESTAssignStuctFlag(StruChkHelper.TESTASSIGNStructFlagSet);


      private final Set<CHECKType> activeChecks;

      CHECKConfig(Set<CHECKType> activeChecks) {
         this.activeChecks = activeChecks;
      }

      public boolean checkIsActive(CHECKType check) {
         return activeChecks.contains(check);
      }
   }

   /** These sets describe supported configurations of
    * StructureChecks.
    */
   private static final Set<CHECKType> ALLCheckSet;
   private static final Set<CHECKType> ALLNonExclusiveSet;
   private static final Set<CHECKType> CHECKStructFlagSet;
   private static final Set<CHECKType> ASSIGNStructFlagSet;

   /** These are for testing only */
   private static final Set<CHECKType> TESTBaseSet;
   private static final Set<CHECKType> TESTCheckStructFlagSet;
   private static final Set<CHECKType> TESTASSIGNStructFlagSet;

   static
   {  ALLCheckSet = EnumSet.allOf(StruChkHelper.CHECKType.class);

      /** set of checks which are not exclusive so that we can add other ones later */
      ALLNonExclusiveSet = EnumSet.copyOf(ALLCheckSet);
      ALLNonExclusiveSet.remove(CHECKType.testComponentNormalizer);
      ALLNonExclusiveSet.remove(CHECKType.assignStructFlag);
      ALLNonExclusiveSet.remove(CHECKType.checkStructFlag);

      /** Add the StructFlag based check */
      CHECKStructFlagSet = EnumSet.copyOf(ALLNonExclusiveSet);
      CHECKStructFlagSet.add(CHECKType.checkStructFlag);

      /** Add check to infer StructFlag */
      ASSIGNStructFlagSet = EnumSet.copyOf(ALLNonExclusiveSet);
      ASSIGNStructFlagSet.add(CHECKType.assignStructFlag);

      TESTBaseSet = EnumSet.copyOf(ALLNonExclusiveSet);
      TESTBaseSet.remove(CHECKType.dbComponentNormalizer);
      TESTBaseSet.add(CHECKType.testComponentNormalizer);

      TESTCheckStructFlagSet = EnumSet.copyOf(TESTBaseSet);
      TESTCheckStructFlagSet.add(CHECKType.checkStructFlag);

      TESTASSIGNStructFlagSet = EnumSet.copyOf(TESTBaseSet);
      TESTASSIGNStructFlagSet.add(CHECKType.assignStructFlag);
   }
}
