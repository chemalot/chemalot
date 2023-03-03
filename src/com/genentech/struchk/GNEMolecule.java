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
package com.genentech.struchk;

import java.util.*;

import com.aestel.utility.Message;
import com.aestel.utility.Message.Level;
import com.genentech.struchk.oeStruchk.OEStruchk;
import com.genentech.struchk.oeStruchk.OEStruchk.StructureFlag;

/**
 * This is a wrapper around {@link OEStruchk} to be used in dataloader.
 * It presets all relevant parameters and provides getter methods to the output.
 * @author albertgo 2008
 *
 */
public class GNEMolecule {
   private final boolean hasError;
   private final int nChiral;
   private final int nChiralSpecified;
   private final int nChiralNonTetrahedral;
   private final int nChiralNonTetrahedralSpecified;
   private final int nNonChiralSp3;  // non chiral sp3 centers with stereo eg C[C@H]1C[C@H]1C
   private final int nNonChiralSp3Specified;
   private final int nStereoDBondSpecified;
   private final int nStereoDBond;
   private final boolean hasChiralFlag;
   private final String inputMf;
   private final String inputMw;
   private final String parentMf;
   private final String parentMw;
   private final String parentMonoIsoMw;
   private final String saltCode;
   private final String saltName;
   private final String saltMW;
   private final String saltMF;
   private final int    saltCount;
   private final String iSmi;
   private final String smi;
   private final String tISmi;
   private final String tSmi;
   private final StructureFlag structFlag;
   private boolean hasStructureFlagError;
   private final List<Message> msgs;
   private final String substanceMolFile;
   private final String substanceISmiles;
   private final String parentMolFile;
   private final String parentStereoMolFile;
   private final String inputMolFile;

   protected GNEMolecule(Builder bldr) {
      this.hasError = bldr.hasError;
      this.inputMolFile = bldr.inputMolFile;
      this.hasChiralFlag = bldr.hasChiralFlag;
      this.saltCode = bldr.saltCode;
      this.saltName = bldr.saltName;
      this.saltMW   = bldr.saltMW;
      this.saltMF   = bldr.saltMF;
      this.saltCount= bldr.saltCount;
      this.inputMw = bldr.inputMw;
      this.inputMf = bldr.inputMf;
      this.parentMw = bldr.parentMw;
      this.parentMonoIsoMw = bldr.parentMonoIsoMw;
      this.parentMf = bldr.parentMf;
      this.substanceMolFile = bldr.substanceMolFile;
      this.substanceISmiles = bldr.substanceISmiles;
      this.parentMolFile = bldr.parentMolFile;
      this.parentStereoMolFile = bldr.parentStereoMolFile;
      this.tISmi = bldr.tISmi;
      this.tSmi = bldr.tSmi;
      this.iSmi = bldr.iSmi;
      this.smi = bldr.smi;
      this.structFlag = bldr.structFlag;
      this.hasStructureFlagError = bldr.hasStructureFLagError;
      this.nChiral = bldr.nChiral;
      this.nChiralSpecified = bldr.nChiralSpecified;
      this.nNonChiralSp3 = bldr.nNonChiralSp3;
      this.nChiralNonTetrahedral = bldr.nChiralNonTetrahedral;
      this.nChiralNonTetrahedralSpecified = bldr.nChiralNonTetrahedralSpecified;
      this.nNonChiralSp3Specified = bldr.nNonChiralSp3Specified;
      this.nStereoDBond = bldr.nStereoDBond;
      this.nStereoDBondSpecified = bldr.nStereoDBondSpecified;
      this.msgs = bldr.msgs;
   }


   /**
    * Return true if there have been errors constructing this molecule.
    *
    * Examples include invalid valences, invalid stereo information.
    */
   public boolean hasError() {
      return hasError;
   }

   /**
    * Get molecular weight of the input structure.
    * Counter ions will not be included.
    */
   public String getInputMolWeight() {
      return inputMw;
   }

   /**
    * Get molecular formula of the input structure.
    * Counter ions will not be included.
    */
   public String getInputMolFormula() {
      return inputMf;
   }

   /**
    * Get molecular formula of the parent structure.
    * Counter ions will not be included.
    */
   public String getParentMolFormula() {
      return parentMf;
   }

   /**
    * Get molecular weight of the parent structure.
    * Counter ions will not be included.
    */
   public String getParentMolWeight() {
      return parentMw;
   }

   /**
    * Get molecular weight of the parent structure computed using the most abundant isotope.
    * Counter ions will not be included.
    */
   public String getParentMonoIsotopicMolWeight() {
      return parentMonoIsoMw;
   }

   /**
    * Get genentech id for the salt of this molecule ("1" for parent = no salt).
    */
   public String getSaltCode() {
      return saltCode;
   }

   /**
    * Get genentech name for salt of this molecule.
    */
   public String getSaltName() {
      return saltName;
   }

   /**
    * Get count of counter ions in input molfile.
    */
   public int getSaltCount() {
      return saltCount;
   }


   /** MW of the salt or "0" if parent */
   public String getSaltMW() {
      return saltMW;
   }

   /** MF of the salt or "" if parent */
   public String getSaltMF() {
      return saltMF;
   }

   /**
    * Return the canonical isomeric smiles of the ParentSubstance.
    * The counter ions have been striped and standardizations applied.
    *
    * @TODO rename to get ParentISmi
    */
   public String getISmi() {
      return iSmi;
   }

   /**
    * Return the canonical non-isomeric smiles of the ParentSubstance.
    * The counter ions have been striped and standardizations applied.
    *
    * @TODO rename to get ParentSmi
    */
   public String getSmi() {
      return smi;
   }

   /**
    * Return canonical smiles of unique tautomer with stereochemistry.
    * This can be used as unique identifier for one parent substance if the
    * substance has the {@link StructureFlag} "No Stereo" or
    * "Single Known Stereoisomer".
    * If the stereochemistry is unknown there might be multiple Substances
    * with the same Structure. These will differ in identity.
    *
    * The tautomer is generated using the heuristic in Quacpac and is not
    * necessarily chemically correct but it is guaranteed to be unique.
    * It should be used for full structure comparisons only and not shown to the user.
    */
   public String getTautomerISmi() {
      return tISmi;
   }

   /**
    * Return canonical smiles of unique tautomer with no stereochemistry.
    * This can be used to look for all isomers of a given compound.
    *
    * The tautomer is generated using the heuristic in Quacpac and is not
    * necessarily chemically correct but it is guaranteed to be unique.
    * It should be used for full structure comparisons only and not shown to the user.
    */
   public String getTautomerSmi() {
      return tSmi;
   }

   /**
    * Return genentech structure flag.
    *
    * This is one of the values returned by {@link StructureFlag#getName}.
    */
   public String getStructureFlag() {
      return structFlag.getName();
   }

   /**
    * If during checking of the structure input flag it was found that the input
    * was inconsistent with the structure drawing this will return true.
    *
    * This is meant to be used with errorAsWarnings.
    */
   public boolean hasStructureFlagError() {
      return hasStructureFlagError;
   }

   /**
    * Return total number of chiral atom with specified or unspecified
    *    r/s stereochemistry including nontetrahedral Chiral centers.
    */
   public int getNChiral() {
      return nChiral;
   }

   /**
    * Return number of chiral atom with specified r/s stereochemistry including nontetrahedral Chiral centers.
    */
   public int getNChiralSpecified() {
      return nChiralSpecified;
   }

   /**
    * Return number of chiral atom that are non tetrahedral eg. atropisomeirc centers.
    */
   public int getNChiralNonTetrahedral() {
      return nChiralNonTetrahedral;
   }


   /**
    * Return number of chiral atom that are non tetrahedral eg. atropisomeirc centers
    * that are specified.
    */
   public int getNChiralNonTetrahedralSpecified() {
      return nChiralNonTetrahedralSpecified;
   }


   /**
    * Return total number of atoms with stereogenic exocyclic bonds which have
    * specified or unspecified stereochemistry.
    *
    * @return eg. 2 for 1,4 dimethylcyclohexan CC1CCC(C)CC1.
    */
   public int getNNonChiralSp3() {
      return nNonChiralSp3;
   }

   /**
    * Return number of atoms with stereogenic exocyclic bonds which have
    * specified stereochemistry.
    *
    * @return eg. 2 for 1,4 trans-dimethylcyclohexan C[C@H]1CC[C@H](C)CC1
    *         and 0 for CC1CCC(C)CC1.
    */
   public int getNNonChiralSp3Specified() {
      return nNonChiralSp3Specified;
   }

   /**
    * Return total number of stereogenic double bonds with and without
    * specified e/z stereochemistry.
    */
   public int getNStereoDBond() {
      return nStereoDBond;
   }

   /**
    * Return number of stereogenic double bonds with specified e/z stereochemistry.
    */
   public int getNStereoDBondSpecified() {
      return nStereoDBondSpecified;
   }

   /**
    * Return true if original molfile had the chiral flag set.
    */
   public boolean hasChiralFlag() {
      return hasChiralFlag;
   }

   /**
    * Return molecule with one copy of counter ion. Standardizations are applied
    * before returning this. Quacpac tautomer normalization will not be applied to this.
    */
   public String getSubstanceMolFile() {
      return substanceMolFile;
   }

   /**
    * Return molecule with one copy of counter ion. Standardizations are applied
    * before returning this. Quacpac tautomer normalization will not be applied to this.
    */
   public String getSubstanceISmiles() {
      return substanceISmiles;
   }

   /**
    * Return molecule from parent molecule normalized with by standardizations.
    * Quacpac tautomer normalization will not be applied to this. Counter Ions
    * will have been removed.
    */
   public String getParentMolFile() {
      return parentMolFile;
   }

   /**
    * Return molecule from parent molecule normalized with by standardizations but
    * including stereochemistry from input.
    * ie. Even Single Unknown Stereochemistry and MI and MD molecules might contain wedges.
    */
   public String getParentStereoMolFile() {
      return parentStereoMolFile;
   }

   /**
    * Return molecule as entered by user, atoms from counter ions will be added.
    */
   public String getInputMolFile() {
      return inputMolFile;
   }

   /**
    * @return a string concatenating all warnings or "" if none.
    */
   public String getWarnings() {
      return getMessagesStr(msgs, EnumSet.of(Level.WARNING));
   }

   /**
    * @return a string concatenating all error messages or "" if none.
    */
   public String getErrors() {
      return getMessagesStr(msgs, EnumSet.of(Level.ERROR));
   }

   /**
    * Return newline separated concatenated list of messages at the specified level.
    * @param levels
    */
   private static String getMessagesStr(List<Message> msgs, EnumSet<Level> levels) {
      StringBuilder sb = new StringBuilder();
      for( Message msg : msgs ) {
         if( levels.contains(msg.getLevel()) )
            sb.append(msg.getText()).append("\n");
      }
      return sb.toString().trim();
   }



   /** Builder for {@link GNEMolecule}s to ease construction while making
    * {@link GNEMolecule} immutable.
    *
    */
   public static class Builder {

      int nChiral;
      int nChiralSpecified;
      int nChiralNonTetrahedral;
      int nChiralNonTetrahedralSpecified;
      int nNonChiralSp3;  // non chiral sp3 centers with stereo eg C[C@H]1C[C@H]1C
      int nNonChiralSp3Specified;
      int nStereoDBondSpecified;
      int nStereoDBond;
      boolean hasChiralFlag;
      String inputMf;
      String inputMw;
      String parentMf;
      String parentMw;
      String parentMonoIsoMw;
      String saltCode;
      String saltName;
      String saltMF;
      String saltMW;
      int    saltCount;
      String iSmi;
      String smi;
      String tISmi;
      String tSmi;
      boolean hasStructureFLagError;
      StructureFlag structFlag;
      List<Message> msgs;
      String parentMolFile;
      String parentStereoMolFile;
      String substanceMolFile;
      String substanceISmiles;
      String inputMolFile;
      boolean hasError;


      GNEMolecule build() {
         return new GNEMolecule(this);
      }

      public Builder setHasError(boolean hasErr) {
         hasError = hasErr;
         return this;
      }

      public Builder setNChiral(int chiral, int chiralSpecified, int nChiralNonTetrahedral, int nChiralNonTetrahedralSpecified) {
         this.nChiral = chiral;
         this.nChiralSpecified = chiralSpecified;
         this.nChiralNonTetrahedral = nChiralNonTetrahedral;
         this.nChiralNonTetrahedralSpecified = nChiralNonTetrahedralSpecified;

         assert nChiral >= nChiralNonTetrahedral :
            "nchiral must include the nonTetraherals chiral centers: " + chiral + " " + nChiralNonTetrahedral;
         assert nChiralSpecified >= nChiralNonTetrahedralSpecified :
            "nchiralSpecified must include the nonTetraherals chiral centers: " + chiralSpecified + " " + nChiralNonTetrahedralSpecified;

         return this;
      }

      public Builder setNNonChiralSp3(int nonChiralSp3, int nonChiralSp3Specified) {
         nNonChiralSp3 = nonChiralSp3;
         nNonChiralSp3Specified = nonChiralSp3Specified;
         return this;
      }
      public Builder setNStereoDBond(int stereoDBond, int stereoDBondSpecified) {
         nStereoDBond = stereoDBond;
         nStereoDBondSpecified = stereoDBondSpecified;
         return this;
      }
      public Builder setHasChiralFlag(boolean hasChiralFlag) {
         this.hasChiralFlag = hasChiralFlag;
         return this;
      }
      public Builder setParentMf(String mf) {
         this.parentMf = mf;
         return this;
      }
      public Builder setParentMw(String mw) {
         this.parentMw = mw;
         return this;
      }
      public Builder setParentMonoIsotpicMw(String mw) {
         this.parentMonoIsoMw = mw;
         return this;
      }
      public Builder setInputMf(String mf) {
         this.inputMf = mf;
         return this;
      }
      public Builder setInputMw(String mw) {
         this.inputMw = mw;
         return this;
      }
      public Builder setSaltCode(String saltCode, String saltName, int saltCount,
                                 String mw, String mf) {
         this.saltCode = saltCode;
         this.saltName = saltName;
         this.saltCount= saltCount;
         this.saltMF   = mf;
         this.saltMW   = mw;
         return this;
      }
      public Builder setISmi(String smi) {
         iSmi = smi;
         return this;
      }
      public Builder setSmi(String smi) {
         this.smi = smi;
         return this;
      }
      public Builder setTSmi(String iSmi, String smi) {
         tSmi = smi;
         tISmi = iSmi;
         return this;
      }
      public Builder setStructFlag(StructureFlag structureFlag) {
         this.structFlag = structureFlag;
         return this;
      }
      public Builder hasStructureFlagError(boolean hasStructureFlagError)
      {  this.hasStructureFLagError = hasStructureFlagError;
         return this;
      }
      public Builder setMsgs(List<Message> msgs) {
         this.msgs = new ArrayList<Message>(msgs);
         return this;
      }
      public Builder setSubstanceMolFile(String substanceMolFile) {
         this.substanceMolFile = substanceMolFile;
         return this;
      }
      public Builder setSubstanceISmiles(String substanceISmiles) {
         this.substanceISmiles = substanceISmiles;
         return this;
      }
      public Builder setParentMolFile(String parentMolFile) {
         this.parentMolFile = parentMolFile;
         return this;
      }
      public Builder setParentStereoMol(String parentStereoMolFile) {
         this.parentStereoMolFile = parentStereoMolFile;
         return this;
      }
      public Builder setInputMolFile(String inputMolFile) {
         this.inputMolFile = inputMolFile;
         return this;
      }

      /** works only after setMsgs was called */
      public void addMsgs(Message msg) {
         this.msgs.add(msg);
      }

      public String getSubstanceMolFile() {
         return substanceMolFile;
      }

      public String getParentMolFile() {
         return parentMolFile;
      }

      /** molfile of parent including original stereochemistry.
       * This even single unknown stereoisomer might contain wedges in this case.
       */
      public String getParentStereoMolFile() {
         return parentStereoMolFile;
      }

      public String getInputMolFile() {
         return inputMolFile;
      }
   }
}
