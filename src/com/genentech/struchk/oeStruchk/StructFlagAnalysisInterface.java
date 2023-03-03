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

import com.genentech.struchk.oeStruchk.OEStruchk.StructureFlag;

/**
 * A StructureCheck class which analyzes the Stereo chemistry and the stereo info
 * and is able to return the correct {@link StructureFlag} flag.
 *
 * @author albertgo
 *
 */
public interface StructFlagAnalysisInterface extends StructureKeeperInterface {

	/** structure Keeper Name for structure Keeper that returns normalized and tautomerized structures
       before setreo centers of non SKS compounds are removed.
       Note: that invalid stereocenters are still removed.
    */
   public static final String STEREONormalizedKeeper = "STEREONormalized";


   /**
    * @return {@link StructureFlag} of last structure passed to {@link #checkStructure}.
    */
   public StructureFlag getStructureFlag();

   /**
    * return true if while checking the flag there was an inconsistency such that
    * the flag is not compatible with the drawing.
    */
   public boolean hasStructureFlagError();

   /**
    * To be called to reset the StereoAnalyser before calling applyRules.
    */
   public void reset();

   /** As of last call to {@link #check} **/
   int getNChiral();

   /** As of last call to {@link #check} **/
   int getNChiralSpecified();

   /** As of last call to {@link #check} **/
   int getNNonChiralStereo();

   /** As of last call to {@link #check} **/
   int getNNonChiralStereoSpecified();

   /** As of last call to {@link #check} **/
   int getNStereoDBond();

   /** As of last call to {@link #check} **/
   int getNStereoDBondSpecified();

   /** return number of chiral centers which are not tetrahedral and were drawn
    * using wedges in last checked molecule.
    *
    * At this time this includes only Atropisomeric centers.
    * As of last call to {@link #check}
    **/
   public int getNChiralNonTetrahedral();

   /** return number of chiral centers which are not tetrahedral and were drawn
    * using wedges in last checked molecule.
    *
    * At this time this includes only Atropisomeric centers.
    * As of last call to {@link #check}
    **/
   public int getNChiralNonTetrahedralSpecified();
}
