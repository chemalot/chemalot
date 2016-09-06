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

import java.io.IOException;
import java.net.URL;
import java.util.HashSet;
import java.util.Set;

import openeye.oechem.OEGraphMol;
import openeye.oechem.oechem;

import org.jdom.JDOMException;

import com.aestel.io.IOUtil;
import com.aestel.utility.DataFormat;
import com.aestel.utility.Message;
import com.aestel.utility.Message.Level;
import com.aestel.utility.exception.IncorrectInputException;
import com.genentech.oechem.tools.OETools;
import com.genentech.struchk.oeStruchk.*;
import com.genentech.struchk.oeStruchk.OEStruchk.StructureFlag;
import com.genentech.struchk.oeStruchk.StruChkHelper.CHECKConfig;
import com.genentech.struchk.oeStruchk.StruChkHelper.CHECKType;

/**
 * This is a wrapper around {@link OEStruchk} to be used in dataloader.
 * It presents all relevant parameters and returns a {@link GNEMolecule} object
 * that provides getter methods to the output.
 *
 * @author albertgo 2008
 *
 */
public class Normalizer {
	public static final String DLCONFIGFile = "substance/ChemRegister.xml";

	/** structure checker which assigns gneStructureFlag */
   private OEStruchk strchkAssignFlag;
   /** structure checker which takes gneStructureFlag and checks validity */
   private OEStruchk strchkCheckFlag;
   private OEGraphMol currentMol = null;

   private int nMessages = 0;
   private int nErrors = 0;
   private int nStruct = 0;

   private static final Set<CHECKType>NOEXCLUSIONS = new HashSet<CHECKType>(0);

   /**
    * Create a Normalizer with default parameters.
    */
   public Normalizer(boolean errorAsWarning)
         throws IncorrectInputException, JDOMException, IOException {
      this(NOEXCLUSIONS, errorAsWarning);
   }
   /**
    * Create a Normalizer with default parameters.
    *
    * Exclusions can be defined like:<br/>
    * <code>Set<CHECKType> exclusions = EnumSet.of(CHECKType.closeAtomCheck);</code>
    *
    * @param exclusions set of rules not to apply while validating.
    */
   public Normalizer(Set<CHECKType> exclusions, boolean errorAsWarning)
            throws IncorrectInputException, JDOMException, IOException {
      URL confFile = IOUtil.getConfigUrl("Struchk.xml", "",
            "/com/genentech/struchk/oeStruchk", false);

      CHECKConfig strChkConfig = CHECKConfig.CHECKStructFlag;
      CHECKConfig strChkConfig2 = CHECKConfig.ASSIGNStructFlag;

      // create OEStruchk from config file
      strchkCheckFlag  = new OEStruchk( confFile, strChkConfig,  exclusions, errorAsWarning );
      strchkAssignFlag = new OEStruchk( confFile, strChkConfig2, exclusions, errorAsWarning );

      currentMol = new OEGraphMol();
   }

   /**
    * Normalize structure from smiles string.
    */
   public GNEMolecule normalizeSmi(String ismi, String gneStructFlag) {
      currentMol.Clear();
      OETools.smiToMol( currentMol, ismi);
      boolean hasChiralFlag   = false;

      OEStruchk checker = strchkAssignFlag;
      StructureFlag sFlag = null;
      String msg = null;
      if(gneStructFlag != null && gneStructFlag.length() > 0) {
         sFlag = StructureFlag.fromString(gneStructFlag);
         if( sFlag != null )
            checker = strchkCheckFlag;
         else
            msg = String.format("Unknown Stereochemistry: %s", gneStructFlag);
      }

      boolean hasError = false;
      if(! checker.applyRules(currentMol, sFlag))
      {  hasError = true;
         nErrors++;
      }

      OEGraphMol outMol = checker.getTransformedMol("parentAllStereo");

      GNEMolecule.Builder gneMolBuilder =
         oeMolToGneMolBuilder(checker, currentMol, outMol, hasChiralFlag, hasError);

      if( msg != null )
         gneMolBuilder.addMsgs( new Message(msg, Level.ERROR, null));

      nStruct++;
      return gneMolBuilder.build();
   }

   /**
    * Strip the mol and return a stripped mol object with associated data in the
    * SD tags.
    *
    * @param molStr input molfile.
    * @param gneStructFlag if null the flag will be assigned assuming unspecified
    *                      centers are mixtures.
    * @return {@link GNEMolecule} containing transformed molecule and associated data.
    */
   public GNEMolecule normalizeMol(String molStr, String gneStructFlag) {
      currentMol.Clear();
      OETools.stringToMol(currentMol, molStr);
      boolean hasChiralFlag   = oechem.OEMDLHasParity( currentMol );

      OEStruchk checker = strchkAssignFlag;
      StructureFlag sFlag = null;
      String msg = null;
      if(gneStructFlag != null && gneStructFlag.length() > 0) {
         sFlag = StructureFlag.fromString(gneStructFlag);
         if( sFlag != null )
            checker = strchkCheckFlag;
         else
            msg  = String.format("Unknown Stereochemistry: %s", gneStructFlag);
      }

      boolean hasError = (msg != null);
      if(! checker.applyRules(molStr, null, sFlag))
      {  hasError = true;
         nErrors++;
      }

      OEGraphMol outMol = checker.getTransformedMol("parentAllStereo");

      GNEMolecule.Builder gneMolBuilder =
               oeMolToGneMolBuilder(checker, currentMol, outMol,hasChiralFlag, hasError)
        .setInputMolFile(molStr);

      if( msg != null )
         gneMolBuilder.addMsgs( new Message(msg, Level.ERROR, null));

      nStruct++;
      return gneMolBuilder.build();
   }


   /**
    * Strip the mol and return a stripped mol object with associated data in the
    * SD tags.
    *
    * @param molStr input molfile.
    * @param gneStructFlag if null the flag will be assigned assuming unspecified
    *                      centers are mixtures.
    * @return {@link GNEMolecule} containing transformed molecule and associated data.
    */
   public GNEMolecule normalizeOEMol(OEGraphMol mol, String gneStructFlag) {
      currentMol.Clear();

      OEStruchk checker = strchkAssignFlag;
      StructureFlag sFlag = null;
      String msg = null;
      if(gneStructFlag != null && gneStructFlag.length() > 0) {
         sFlag = StructureFlag.fromString(gneStructFlag);
         if( sFlag != null )
            checker = strchkCheckFlag;
         else
            msg = String.format("Unknown Stereochemistry: %s", gneStructFlag);
      }

      boolean hasError = false;
      if(! checker.applyRules(mol, sFlag))
      {  hasError = true;
         nErrors++;
      }

      OEGraphMol outMol = checker.getTransformedMol("parentAllStereo");

      GNEMolecule.Builder gneMolBuilder =
               oeMolToGneMolBuilder(checker, currentMol, outMol, false, hasError)
        .setInputMolFile(OETools.molToString(mol));

      if( msg != null )
         gneMolBuilder.addMsgs( new Message(msg, Level.ERROR, null));

      nStruct++;
      return gneMolBuilder.build();
   }


   private GNEMolecule.Builder oeMolToGneMolBuilder(OEStruchk checker,
         OEGraphMol inMol, OEGraphMol parntMol, boolean hasChiralFlag, boolean hasError) {
      GNEMolecule.Builder gneMolBuilder = new GNEMolecule.Builder()
         .setHasError(hasError)
         .setHasChiralFlag(hasChiralFlag)
         .setInputMf( oechem.OEMolecularFormula(inMol))
         .setInputMw( DataFormat.formatNumber(oechem.OECalculateMolecularWeight(inMol,true), "r2"))
         .setParentMf(oechem.OEMolecularFormula(parntMol))
         .setParentMw(DataFormat.formatNumber(oechem.OECalculateMolecularWeight(parntMol,true), "r2"))
         .setParentMonoIsotpicMw(DataFormat.formatNumber(OETools.getIsotpictWeight(parntMol), "r4"))
         .setSaltCode(        checker.getSaltCode(), checker.getSaltName(), checker.getSaltCount(),
                              checker.getSaltMW(), checker.getSaltMF())
         .setTSmi(            checker.getTransformedIsoSmiles("tautomer"),
                              checker.getTransformedSmiles("tautomer"))
         .setISmi(            checker.getTransformedIsoSmiles("parent"))
         .setSmi(             checker.getTransformedSmiles("parent"))
         .setStructFlag(      checker.getStructureFlag())
         .hasStructureFlagError( checker.hasStructureFlagError() )
         .setMsgs(            checker.getStructureMessages(null))
         .setParentStereoMol( checker.getTransformedMolfile("parentAllStereo"))
         .setParentMolFile(   checker.getTransformedMolfile("parent"))
         .setSubstanceMolFile(checker.getTransformedMolfile("substance"))
         .setSubstanceISmiles(checker.getTransformedIsoSmiles("substance"))
         .setNChiral(checker.getNChiral(), checker.getNChiralSpecified())
         .setNNonChiralSp3(checker.getNNonChiralSp3(), checker.getNNonChiralSp3Specified())
         .setNStereoDBond(checker.getNStereoDBond(), checker.getNStereoDBondSpecified());

      return gneMolBuilder;
   }


   /** after calling delete this Normalizer should not be used again */
   public void close() {
      currentMol.delete();

      strchkAssignFlag.delete();
      strchkCheckFlag.delete();
   }

   int getNMessages() {
      return nMessages;
   }

   int getNErrors() {
      return nErrors;
   }

   int getNStruct() {
      return nStruct;
   }

   public static void main(String ... args) throws IncorrectInputException, JDOMException, IOException {
      Normalizer normalizer = new Normalizer(false);
      String molStr;
      molStr = IOUtil.fileToString("c:/tmp", "t.sdf");
//      GNEMolecule gMol = normalizer.normalizeSmi("CCOC.Cl", "No Stereo");
      GNEMolecule gMol = normalizer.normalizeMol(molStr, "Single Known Stereoisomer");

      if(gMol.hasError())
         System.err.println(gMol.getErrors());

      // even if errors occur oeStructChk will try to normalize the structure
      System.err.printf("sFlag=%s\n",   gMol.getStructureFlag());

      // to be used for tautomeric uniqueness check including stereo chemistry
      System.err.printf("tiSmi =%s\n",   gMol.getTautomerISmi());

      // to be used for tautomeric uniqueness check excluding stereo chemistry
      System.err.printf("tSmi =%s\n",   gMol.getTautomerSmi());

      // to be used for display of the parent structure
      System.err.printf("iSmi =%s\n",   gMol.getParentMolFile());

      System.err.printf("salt =%s * %s\n", gMol.getSaltCount(), gMol.getSaltCode());

      // To be used for display of the substance: (only one copy of counter ion)
      System.err.printf("SubstISmi=%s\n", gMol.getSubstanceMolFile());

      // free up resources
      normalizer.close();
   }
}
