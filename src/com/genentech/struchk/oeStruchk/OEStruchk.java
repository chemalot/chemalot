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

import java.io.*;
import java.net.URL;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import openeye.oechem.*;

import org.apache.commons.cli.*;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.JDOMException;
import org.jdom.input.SAXBuilder;

import com.aestel.utility.Message;
import com.aestel.utility.MessageList;
import com.aestel.utility.WarnOnlyMessageList;
import com.genentech.oechem.tools.AtomHasPropertyFunctor;
import com.genentech.oechem.tools.OETools;
import com.genentech.struchk.oeStruchk.StruChkHelper.CHECKConfig;
import com.genentech.struchk.oeStruchk.StruChkHelper.CHECKType;
import com.genentech.struchk.oeStruchk.StructureCheckInterface.HydrogenMode;

/**
 * Business rule checker for chemical structures.
 *
 * This implementation is based on the open eye toolkit and uses an xml
 * configuration file to define the transformations and unwanted structures.
 *
 */
public class OEStruchk {

   private static final boolean PRINTOEObjects = false;

   public static enum StructureFlag {
      NOStereo("No Stereo", "NS"),
      SINGLEStereoisomer("Single Known Stereoisomer", "Single Stereoisomer", "SI"),
      SINGLEUnknownStereoIsomer("Single Unknown Stereoisomer", "SUI"),
      MIXTUREOfEnantiomers("Mixture of Enantiomers", "ME"),
      MIXTUREOfDiastereomers("Mixture of Diastereomers", "Mixture of Stereoisomers", "MI"),
      UNCERTAINStructure("Uncertain Structure");

      private final String[] names;
      private static final Map<String, StructureFlag> nameToSF
                                          = new HashMap<String, StructureFlag>();
      static {
         for(StructureFlag si : StructureFlag.values())
            for(String name : si.names)
               nameToSF.put(name, si);
      }

      private StructureFlag(String... names) {
         this.names = names;
      }

      public String getName() {
         return names[0];
      }

      public static StructureFlag fromString(String name) {
         StructureFlag sf = nameToSF.get(name);
         return sf;
      }
   }


   /** tag to mark atoms in oechem according to their index */
   public static final int ATOMNumTag = oechem.OEGetTag("STRCHKANum");

   /** tag to store atomic coordinates of explicit hydrogens */
   static final int ATOMExpHX = oechem.OEGetTag("STRCHKAExpHX");
   static final int ATOMExpHY = oechem.OEGetTag("STRCHKAExpHY");
   static final int ATOMExpHZ = oechem.OEGetTag("STRCHKAExpHZ");


   /** tag to mark bonds and atoms on which the stereochemistry was removed
    * by {@link DoubleBondStereoClean}. */
   public static final int STEREOClearTag = oechem.OEGetTag("STRCHStereoClear");

   /** tag to mark Bridgeheads which can not be inverted independently as in
    *  C[C@H]1C[C@@H]2[C@@H](C2)C1 compare for doc-files\BycycleStereo.pdf. */
   public static final int ISStrainedBridgeHead = oechem.OEGetTag("STRCHIsStrainedBridgeHead");

   /** tag atoms as chiral which are not correctly recognized as such by OEChem
    * eg.  F[C@@H]([3H])C */
   public static final int ISChiralNotRecognized = oechem.OEGetTag("STRCHIsChiralNotRecognized");

   /** tag to mark bonds and atoms on had stereochemistry in the before tautomerization. */
   public static final int HADStereoTag = oechem.OEGetTag("HADStereoTag");

   /** tag to mark sp3 atoms which are in rings and therefore could bear stereo information
    * eg. in CC1CCC(C)CC1. Used by {@link CheckStructureFlag}. */
   public static final int NONChiralStereoAtomTag  = oechem.OEGetTag("NONChiralStereo");
   /** flag on whole molecule flagging a previous call to {@see #flagNonchiralStereoAtoms}*/
   public static final int NONChiralStereoAssignedTag = oechem.OEGetTag("NONChiralStereoMol");

   /** flag on atoms and bonds if they are detected to be part of an atropisomeric cneter */
   public static final int ATROPIsomericCenter = oechem.OEGetTag("AtropIsomeric");



   private static final Pattern thickBondPattern = Pattern.compile("Display=[\"']Bold[\"']");
   //private static final Pattern wigglyBondPattern = Pattern.compile("[\\n\\r]+.+[\\n\\r]+.+[\\n\\r]+.+[\\n\\r]+.+[\\n\\r]+[ \\d]{6}  1  4[ \\n\\r]");
   private static final Pattern atomLabelPattern = Pattern.compile("M  STY.*SUP");

   /** Set of rules to be applied */
   private final StructureCheckInterface[] rules;

   /** current molecule */
   private final OEGraphMol oeMol = new OEGraphMol();

   /** temporary molecule object to be used within one method at the time to
    * save garbage collection cost */
   private final OEGraphMol tmpMol = new OEGraphMol();
   private final MyCannonicalize canonizer = new MyCannonicalize();

   /** Functors create only once to save garbage collection */
   private OEHasAtomicNum hAtomFunct = new OEHasAtomicNum(1); // find explicit H
   private AtomHasPropertyFunctor atomHadExplicitHFunct = new AtomHasPropertyFunctor(ATOMExpHX);

   private final MessageList structMessages;

   private final String checkForThickBondDesc;
   //private final String checkForWigglyBondDesc;
   private final String checkForAtomLabelDesc;

   private final ComponentNormalizerInterface componentNormalizer;
   private final StructFlagAnalysisInterface sflagChecker;

   private final Map<String,StructureKeeperInterface> keeperMap
                                = new HashMap<String,StructureKeeperInterface>();

   /** stores the chiral flag of the current molecule */
   private boolean hasChiralFlag;

   /** true if the current molecule is empty and was replaced by "*" */
   private boolean isEmptyMol;

   private static final Set<CHECKType>NOEXCLUSIONS = new HashSet<CHECKType>(0);


   /** Create a structure checking engine from the xml configuration file
    *
    * @throws IOException on reading config file
    * @throws JDOMException on problems in the xml file */
   public OEStruchk(URL cFile, CHECKConfig config, boolean errorsAreWarnings) throws JDOMException, IOException {
      this(cFile, config, NOEXCLUSIONS, errorsAreWarnings);
   }
   /** Create a structure checking engine from the xml configuration file
    * excluding a subset of the rules.
    *
    * Exclusions can be defined like:<br/>
    * <code>Set<CHECKType> exclusions = EnumSet.of(CHECKType.closeAtomCheck);</code>
    *
    * @param exclusions set of rules to be excluded from execution.
    *
    * @throws IOException on reading config file
    * @throws JDOMException on problems in the xml file */
   public OEStruchk(URL cFile, CHECKConfig config, Set<CHECKType> exclusions, boolean errorsAreWarnings)
                   throws JDOMException, IOException {

      SAXBuilder builder = new SAXBuilder();
      Document confFile = builder.build(cFile);

      String reqOEChemVersion=confFile.getRootElement().getAttributeValue("oechemVersion");
      if(! reqOEChemVersion.contains(Integer.toString(oechem.OEChemGetVersion())))
         throw new Error(String.format("Wrong OEchem Version: required=%s actual=%d\n"
                              +"Change Struchk.xml or install correct version\n",
                     reqOEChemVersion, oechem.OEChemGetVersion()));

      if( ! errorsAreWarnings )
         structMessages = new MessageList();
      else
         structMessages = new WarnOnlyMessageList();

      List<StructureCheckInterface> ruleList = new ArrayList<StructureCheckInterface>();

      String _checkForThickBonds = null;
      //String _checkForWigglyBonds= null;
      String _checkForAtomLabels  = null;
      ComponentNormalizerInterface normalizer = null;
      StructFlagAnalysisInterface sflagChk = null;
      FlagNonChiralStereoCenters stereoFlagger = null;

      HydrogenMode hMode = HydrogenMode.ANY;
      // read list of rules from xml file
      // factory method might be better
      printOEObjects("before rules");
      for(Object ruleElementO : confFile.getRootElement().getChildren()) {
         Element ruleElement = (Element)ruleElementO;

         CHECKType checkType = CHECKType.valueOf(ruleElement.getAttributeValue("id"));
         if(! config.checkIsActive(checkType)) continue;
         if( exclusions.contains( checkType )) continue;

         StructureCheckInterface rule = null;

         if("transform".equals(ruleElement.getName())) {
            rule = new Transformer(ruleElement);

         } else if("badSubstructure".equals(ruleElement.getName())) {
            rule = new BadSubstructureCheck(ruleElement);

         } else if("checkChiral".equals(ruleElement.getName())) {
            rule = new ChiralityCheck(ruleElement);

         } else if("checkDoubleBond".equals(ruleElement.getName())) {
            rule = new DoubleBondCheck(ruleElement);

         } else if("clearBondStereo".equals(ruleElement.getName())) {
            rule = new DoubleBondStereoClean(ruleElement);

         } else if("checkAtomtypes".equals(ruleElement.getName())) {
            rule = new AtomTypeCheck(this, ruleElement);

         } else if("valenceCheck".equals(ruleElement.getName())) {
            rule = new AtomValenceCheck(ruleElement);

         } else if("closeAtomsCheck".equals(ruleElement.getName())) {
            rule = new CloseAtomsCheck(ruleElement);

         } else if("flagNonChiralAtoms".equals(ruleElement.getName())) {
            stereoFlagger = new FlagNonChiralStereoCenters(ruleElement);
            rule = stereoFlagger;

         } else if("twoDCheck".equals(ruleElement.getName())) {
            rule = new TwoDCheck(ruleElement);

         } else if("cleanReactionCenter".equals(ruleElement.getName())) {
            rule = new CleanReactionCenter(ruleElement);

         } else if("tautomerize".equals(ruleElement.getName())) {
            rule = new TautomerStandardizer(ruleElement);

         } else if("wigglyBondCheck".equals(ruleElement.getName())) {
            rule = new WigglyBondCheck(ruleElement);

         } else if("removeHydrogens".equals(ruleElement.getName())) {
            rule = new HydrogenRemover(ruleElement);
            hMode = HydrogenMode.SUPRRESSED;

         } else if("assignStructFlag".equals(ruleElement.getName())) {
            sflagChk = new AssignStructureFlag(ruleElement, stereoFlagger);
            rule = sflagChk;

         } else if("checkStructFlag".equals(ruleElement.getName())) {
            sflagChk = new CheckStructureFlag(ruleElement, stereoFlagger);
            rule = sflagChk;

         } else if("componentNormalizer".equals(ruleElement.getName())) {
            if( normalizer != null )
               throw new Error("Only one ComponentNormalizer may be defined!");
            ComponentNormalizer cNormalizer = new ComponentNormalizer(ruleElement);
            if( cNormalizer.getKeeperName() != null)
               keeperMap.put(cNormalizer.getKeeperName(), cNormalizer);
            normalizer = cNormalizer;
            rule       = cNormalizer;

         } else if("atomLabelCheck".equals(ruleElement.getName())) {
            _checkForAtomLabels = ruleElement.getChildTextTrim("description");
            assert _checkForAtomLabels != null : "Description expected";
            rule = null;

         } else if("thickBondCheck".equals(ruleElement.getName())) {
            _checkForThickBonds = ruleElement.getChildTextTrim("description");
            assert _checkForThickBonds != null : "Description expected";
            rule = null;

//         } else if("wigglyBondCheck".equals(ruleElement.getName())) {
//            _checkForWigglyBonds = ruleElement.getChildTextTrim("description");
//            assert _checkForWigglyBonds != null : "Description expected";
//            rule = null;

         } else if("keepStructure".equals(ruleElement.getName())) {
            rule = new StructureKeeper(ruleElement);
            keeperMap.put(((StructureKeeper)rule).getKeeperName(), (StructureKeeper)rule);

         } else
            throw new Error("Unknown rule element:" + ruleElement.getName());

         if(rule == null ) continue;

         if(rule.getRequiredHydrogenMode() != HydrogenMode.ANY
            && rule.getRequiredHydrogenMode() != hMode )
            throw new Error(
                  String.format("Rule needs %s hydrogen mode but current mode is %s.",
                        rule.getRequiredHydrogenMode().toString(), hMode));

         ruleList.add(rule);
         printOEObjects(rule.getCheckName());
      }
      checkForThickBondDesc = _checkForThickBonds;
      checkForAtomLabelDesc = _checkForAtomLabels;
      componentNormalizer  = normalizer;
      sflagChecker         = sflagChk;

      rules = ruleList.toArray(new StructureCheckInterface[ruleList.size()]);
   }

   private static final Pattern EMPTYMolPattern = Pattern.compile(
            "^.*(\\n\\r*|\\r\\n*).*(\\n\\r*|\\r\\n*).*(\\n\\r*|\\r\\n*)  0  0");

   /**
    * Apply all roles to the input molecule, append any messages to the msgs list.
    *
    * @param mol input molecule which will be modified.
    * @param cdxml input molecule as cdxml, to be used for thick bond check, may be null.
    * @param structFlag StructureFlag for mol, null if mode is EXTERNAL
    * @return false if any errors occur.
    */
   public boolean applyRules(String mol, String cdxml, StructureFlag structFlag ) {
      boolean success = true;

      resetState();

      for( StructureKeeperInterface sKeep : keeperMap.values())
         sKeep.delete();

      if( EMPTYMolPattern.matcher(mol).find()) {
// Disabled empty mol check on request from Gina with OK from ROhan and Fred 2010/02/12
//         structMessages.addMessage(new Message(
//               "Empty Molecules are not allowed.", Message.Level.ERROR, null));
//         success = false;
         OETools.smiToMol(oeMol, "*");
         isEmptyMol = true;

      }else {
         OETools.stringToMol(oeMol, mol);
      }

      hasChiralFlag   = oechem.OEMDLHasParity( oeMol );

      // check for atom labels in the molfile
      if(checkForAtomLabelDesc != null && atomLabelPattern.matcher(mol).find()) {
         structMessages.addMessage(new Message(
               "Functional group abbreviations are not allowed.", Message.Level.ERROR, null));
         success = false;
      }

      // check for thick bonds in the chemdraw xml file
      if(checkForThickBondDesc != null && cdxml != null
            && thickBondPattern.matcher(cdxml).find() ) {
         structMessages.addMessage(new Message(
               "Thick bonds are not allowed.", Message.Level.ERROR, null));
         success = false;
      }

      // check for wiggly bonds in the mol file, oechem ignores them
//      if(checkForWigglyBondDesc != null
//            && wigglyBondPattern.matcher(mol).find() ) {
//         structMessages.addMessage(new Message(
//               "Wiggly bonds are not allowed.", Message.Level.ERROR, null));
//         success = false;
//      }

      return applyRules(success, structFlag );
   }

   /**
    * Apply all roles to the input molecule, append any messages to the msgs list.
    *
    * @param mol input molecule which will not be modified.
    * @return false if any errors occur.
    */
   public boolean applyRules(OEGraphMol mol, StructureFlag inSFlag)
   {  boolean success = true;

      resetState();

      hasChiralFlag   = oechem.OEMDLHasParity( mol );

      for( StructureKeeperInterface sKeep : keeperMap.values())
         sKeep.delete();

      oechem.OEAddMols(oeMol, mol);
      oeMol.SetTitle( mol.GetTitle() );

      return applyRules(success, inSFlag);
   }

   /**
    * Apply all rules to internally stored oeMol.
    * @param structFlag StructureFlag for this molecule, may be null for StruCheckMode = EXTERNAL.
    *
    * @return false on error, all keepers will be set to empty molecules if
    *         a fatal error happens.
    */
   private boolean applyRules(boolean success, StructureFlag structFlag) {
      if(! oeMol.IsValid()) {
         structMessages.addMessage(new Message("Invalid molecule.",Message.Level.ERROR, null));

         // add empty molecules to all keepers
         oeMol.Clear();
         for( StructureKeeperInterface sKeep : keeperMap.values())
            sKeep.checkStructure(oeMol, StructureFlag.NOStereo, structMessages);
         return false;
      }

      saveExplicitHydrogenLocations();

      // tag atoms so that we can distinguish them afterwards
      initMolForChecks(oeMol);
      oechem.OEPerceiveChiral(oeMol);

      // apply all rules to molecule
      for(StructureCheckInterface rule : rules) {
//System.err.printf("rule %s mol %s\n",rule.getDescription(),oechem.OECreateCanSmiString(oeMol));
         success = success & rule.checkStructure(oeMol, structFlag, structMessages);
//System.err.printf("%s %s\n", rule.toString(), OETools.molToCanSmi(oeMol, true));
      }
//System.err.println(OETools.molToString(oeMol));

      oechem.OESuppressHydrogens(oeMol,false,false,true);
      oechem.OEPerceiveChiral(oeMol);

      return success;
   }

   /**
    * Save the 3D location in of explicit Hydrogens so that they can be recreated
    * in molfiles for better stereo representation.
    */
   private void saveExplicitHydrogenLocations() {
      OEAtomBaseIter aIt = oeMol.GetAtoms(hAtomFunct);
      double[] coor = new double[3];
      while( aIt.hasNext() ) {
         OEAtomBase hAt = aIt.next();
         if( hAt.GetExplicitDegree() != 1 ) continue; // more than one bond?

         OEBondBaseIter bIt = hAt.GetBonds();
         OEBondBase bd = bIt.next();
         OEAtomBase neigh = bd.GetNbr(hAt);
         bIt.delete();

         if( neigh.GetImplicitHCount() > 0 || neigh.GetExplicitHCount() > 1 )
            continue; // we are only interested in chiral H, must be unique

         oeMol.GetCoords(hAt, coor);
         neigh.SetDoubleData(ATOMExpHX, coor[0]);
         neigh.SetDoubleData(ATOMExpHY, coor[1]);
         neigh.SetDoubleData(ATOMExpHZ, coor[2]);
         neigh.SetIntData(OEProperty.BondStereo, bd.GetIntData(OEProperty.BondStereo));
      }
      aIt.delete();
   }

   /**
    * Return the transformed molfile of the molecule which was last passed to
    *    {@link #applyRules} at the stage of the named structure keeper.
    *
    * This method will try to add explicit hydrogens with original coordinates
    * to atoms which had explicit hydrogens on input.
    */
   public String getTransformedMolfile(String keeperName) {
      tmpMol.Clear();
      oechem.OEAddMols(tmpMol, getTransformedMol(keeperName));
      OEAtomBaseIter aIt = tmpMol.GetAtoms( atomHadExplicitHFunct );
      double[] coor = new double[3];
      while( aIt.hasNext() ) {
         OEAtomBase at = aIt.next();
         if( at.GetImplicitHCount() != 1 ) continue;

         coor[0] = at.GetDoubleData(ATOMExpHX);
         coor[1] = at.GetDoubleData(ATOMExpHY);
         coor[2] = at.GetDoubleData(ATOMExpHZ);
         int upDownType = at.GetIntData(OEProperty.BondStereo);

         OEAtomBase hAt = tmpMol.NewAtom(1);
         tmpMol.SetCoords(hAt, coor);
         at.SetImplicitHCount(0);

         OEBondBase bd = tmpMol.NewBond(at, hAt, 1);
         bd.SetIntData(OEProperty.BondStereo, upDownType);
      }
      aIt.delete();

      oechem.OEMDLPerceiveParity(tmpMol);
      oechem.OEMDLCorrectBondStereo(tmpMol);

      return OETools.molToString(tmpMol);
   }

   /**
    * return the transformed OEGraphmol at the stage of the structure keeper
    * named keeperName.
    *
    * @param keeperName the name of one of the keepStructure elements in the
    * configuration file, if keeperName is null the final molecule after all
    * transformations is returned.
    *
    *  @return OEGraphMol of the molecule at the stage named keeperName of this
    *    OEStruchk run, NOTE this object will be deleted during the next call to
    *    {@link #applyRules}.
    */
   public OEGraphMol getTransformedMol(String keeperName) {
      if( keeperName == null ) return oeMol;

      assert keeperMap.containsKey(keeperName) : "Unknown keeper " + keeperName;
      return keeperMap.get(keeperName).getMolecule();
   }

   /**
    * Return the transformed canonical isomeric smiles of the molecule which was
    * last passed to {@link #applyRules} at the stage of the structure
    * keeper named keeperName.
    *
    * @param keeperName the name of one of the keepStructure elements in the
    * configuration file, if keeperName is null the final molecule after all
    * transformations is returned.
    *
    * @return may return null if input structure was completely invalid.
    */
   public String getTransformedIsoSmiles(String keeperName) {
      return canonizer.myCanonicalSmi(getTransformedMol(keeperName), true);
   }

   /**
    * Return the transformed canonical smiles of the molecule which was
    * last passed to {@link #applyRules} at the stage of the structure
    * keeper named keeperName.
    *
    * @param keeperName the name of one of the keepStructure elements in the
    * configuration file, if keeperName is null the final molecule after all
    * transformations is returned.
    */
   public String getTransformedSmiles(String keeperName) {
      // we could use MyCannonicalize.myCanonicalSmi(getTransformedMol(keeperName), false)
      // but this is a lot slower, in order to prevent some very esoteric problems
      // with the OE canonicalization algorithm. Since this is less important
      // for the non-isomeric smiles we use the faster direct method.
      return OETools.molToCanSmi(getTransformedMol(keeperName), false);
   }

   /** return number of chiral atoms in last checked molecule */
   public int getNChiral()
   {  return sflagChecker.getNChiral();
   }

   /** return number of chiral atoms wiht specified stereo in last checked molecule */
   public int getNChiralSpecified()
   {  return sflagChecker.getNChiralSpecified();
   }

   /** return number of non-chiral sp3 atoms in last checked molecule */
   public int getNNonChiralSp3()
   {  return sflagChecker.getNNonChiralStereo();
   }

   /** return number of non-chiral sp3 atoms with specified stereo in last checked molecule */
   public int getNNonChiralSp3Specified()
   {  return sflagChecker.getNNonChiralStereoSpecified();
   }

   /** return number of stereo double bonds in last checked molecule */
   public int getNStereoDBond()
   {  return sflagChecker.getNStereoDBond();
   }

   /** return number of  stereo double bonds with specified stereo in last checked molecule */
   public int getNStereoDBondSpecified()
   {  return sflagChecker.getNStereoDBondSpecified();
   }


   /**
    * Get the instance of the {@link StructFlagAnalysisInterface} used to check and
    * analyze the input structure.
    */
   public StructFlagAnalysisInterface getSFlagAnalyzer() {
      return sflagChecker;
   }

   /**
    * @return true if the last molecule passed to {@link #applyRules} had the MDL chiral flag set.
    */
   public boolean hasChiralFlag() {
      return hasChiralFlag;
   }

   /**
    * return Messages generated during the normalization of the last
    * molecule up to the stage of the structure keeper named keeperName.
    *
    * @param keeperName the name of one of the keepStructure elements in the
    * configuration file, if keeperName is null the final list of messages
    * after all transformations is returned.
    *
    * @return a (possibly empty (length == 0)) list of messages.
    */
   public List<Message> getStructureMessages(String keeperName) {
      if( keeperName == null ) return structMessages.getMessages();

      return keeperMap.get(keeperName).getStructureMessages().getMessages();
   }

   /**
    * Return the salt code of the molecule which was
    * last passed to {@link #applyRules} .
    *
    */
   public String getSaltCode() {
      assert componentNormalizer != null : "Componentnormalizer Check not configured!";

      String saltCode = componentNormalizer.getSaltCode();

      return saltCode;
   }

   /**
    * Return the salt name of the molecule which was
    * last passed to {@link #applyRules} .
    *
    */
   public String getSaltName() {
      assert componentNormalizer != null : "Componentnormalizer Check not configured!";

      String saltCode = componentNormalizer.getSaltName();

      return saltCode;
   }

   /** Return the number of repetitions of the counter ion (salt) present in the input
    * molecule.
    */
   public int getSaltCount() {
      assert componentNormalizer != null : "Componentnormalizer Check not configured!";

      return componentNormalizer.getSaltCount();
   }

   /** Return the molecular weight of the counter ion (salt) present in the input
    * molecule or "0" if none.
    */
   public String getSaltMW() {
      assert componentNormalizer != null : "Componentnormalizer Check not configured!";

      return componentNormalizer.getSaltMW();
   }

   /** Return the molecular formula of the counter ion (salt) present in the input
    * molecule or "" if none.
    */
   public String getSaltMF() {
      assert componentNormalizer != null : "Componentnormalizer Check not configured!";

      return componentNormalizer.getSaltMF();
   }

   private void resetState() {
      oeMol.Clear();
      structMessages.clear();
      if(componentNormalizer != null) componentNormalizer.reset();
      if(sflagChecker       != null) sflagChecker.reset();
      isEmptyMol = false;
   }

   public static void main3(String[] args) throws JDOMException, IOException {
      URL url = getResourceURL(OEStruchk.class,"Struchk.xml");

      OEGraphMol mol = new OEGraphMol();
      oechem.OEParseSmiles(mol, "N/C=N/[H]");

      // create OEStruchk from config file
      OEStruchk checker =  new OEStruchk(url, CHECKConfig.TESTAssignStuctFlag, false);
      if(! checker.applyRules(mol, null))
         System.err.println("Unsuccesful");

      System.err.println(checker.getTransformedIsoSmiles(null)
            + " " + checker.getSaltCode() + " " + checker.getStructureFlag());

      List<Message> msgs = checker.getStructureMessages(null);
      for(Message msg : msgs) {
         System.err.println(msg.getText());
      }

   }


   /**
    * Command line interface to {@link OEStruchk}.
    */
   public static void main(String [] args) throws ParseException, JDOMException, IOException {
      long start = System.currentTimeMillis();
      int nMessages = 0;
      int nErrors = 0;
      int nStruct = 0;
      System.err.printf("OEChem Version: %s\n", oechem.OEChemGetVersion());

      // create command line Options object
      Options options = new Options();
      Option opt = new Option("f", true, "specify the configuration file name");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("noMsg", false, "Do not add any additional sd-tags to the sdf file");
      options.addOption(opt);

      opt = new Option("printRules", true, "Print HTML listing all the rules to filename.");
      options.addOption(opt);

      opt = new Option("errorsAsWarnings", false, "Treat errors as warnings.");
      options.addOption(opt);

      opt = new Option("stopForDebug", false, "Stop and read from stdin for user tu start debugger.");
      options.addOption(opt);


      CommandLineParser parser = new PosixParser();
      CommandLine cmd = parser.parse( options, args);
      args = cmd.getArgs();

      if( cmd.hasOption("stopForDebug"))
      {  BufferedReader localRdr = new BufferedReader(new InputStreamReader( System.in ));
         System.err.print("Please press return:");

         localRdr.readLine();
      }

      URL confFile;
      if( cmd.hasOption( "f") )
      {  confFile = new File(cmd.getOptionValue("f")).toURI().toURL();
      } else
      {  confFile = getResourceURL(OEStruchk.class,"Struchk.xml");
      }
      boolean errorsAsWarnings = cmd.hasOption("errorsAsWarnings");

      if(cmd.hasOption("printRules")) {
         String fName = cmd.getOptionValue("printRules");
         PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(fName)));
         OEStruchk structFlagAssigner =  new OEStruchk(confFile,
                                          CHECKConfig.ASSIGNStructFlag, errorsAsWarnings);
         structFlagAssigner.printRules(out);
         out.close();
         return;
      }

      if(args.length < 1) {
         HelpFormatter formatter = new HelpFormatter();
         formatter.printHelp( "oeStruck", options );
         throw new Error( "missing input file\n");
      }

      BufferedReader in = null;
      try {
         in = new BufferedReader(new FileReader(args[0]));
         StringBuilder mol = new StringBuilder();
         StringBuilder data = new StringBuilder();
         PrintStream out = System.out;
         if(args.length > 1) out = new PrintStream(args[1]);

         // create OEStruchk from config file
         OEStruchk structFlagAssigner =  new OEStruchk(confFile,
                                                     CHECKConfig.ASSIGNStructFlag, errorsAsWarnings);

         OEStruchk structFlagChecker = new OEStruchk(confFile,
                                                     CHECKConfig.CHECKStructFlag, errorsAsWarnings);

         Pattern sFlagPat = Pattern.compile("<StructFlag>\\s*([^\\n\\r]+)");
         String line;
         boolean inMolFile = true;
         boolean atEnd = false;
         while(! atEnd ) {
            if((line=in.readLine()) == null) {
               if("".equals(mol.toString().trim())) break;

               if(! inMolFile)
                  throw new Error("Invalid end of sd file!");

               line = "$$$$";
               atEnd = true;
            }

            if( line.startsWith("$$$$")) {

               OEStruchk oeStruchk;
               StructureFlag sFlag = null;
               Matcher mat = sFlagPat.matcher(data);
               if(! mat.find()) {
                  oeStruchk = structFlagAssigner;
               }else {
                  oeStruchk = structFlagChecker;
                  sFlag = StructureFlag.fromString(mat.group(1));
               }
               if(! oeStruchk.applyRules(mol.toString(), null, sFlag))
                  nErrors++;

               out.print(oeStruchk.getTransformedMolfile(null));

               out.print(data);

               if(! cmd.hasOption("noMsg")) {
                  List<Message> msgs = oeStruchk.getStructureMessages(null);
                  if(msgs.size() > 0) {
                     nMessages += msgs.size();

                     out.println("> <errors_oe2>");
                     for( Message msg : msgs )
                        out.printf("%s: %s\n", msg.getLevel(), msg.getText());
                     out.println();
                  }
//System.err.println(oeStruchk.getTransformedMolfile("substance"));
                  out.printf("> <outStereo>\n%s\n\n", oeStruchk.getStructureFlag().getName());
                  out.printf("> <TISM>\n%s\n\n", oeStruchk.getTransformedIsoSmiles(null));
                  out.printf("> <TSMI>\n%s\n\n", oeStruchk.getTransformedSmiles(null));
                  out.printf("> <pISM>\n%s\n\n", oeStruchk.getTransformedIsoSmiles("parent"));
                  out.printf("> <salt>\n%s\n\n", oeStruchk.getSaltCode());
                  out.printf("> <stereoCounts>\n%s.%s\n\n",
                        oeStruchk.countChiralCentersStr(),
                        oeStruchk.countStereoDBondStr());
               }

               out.println(line);
               nStruct++;

               mol.setLength(0);
               data.setLength(0);

               inMolFile = true;
            }else if( ! inMolFile || line.startsWith(">")){
               inMolFile = false;
               data.append(line).append("\n");
            }else {
               mol.append(line).append("\n");
            }
         }

         structFlagAssigner.delete();
         structFlagChecker.delete();

      } catch (Exception e) {
         throw new Error(e);
      }finally {
         System.err.printf("Checked %d structures %d errors, %d messages in %dsec\n",
               nStruct, nErrors, nMessages, (System.currentTimeMillis()-start)/1000);
         if( in != null ) in.close();
      }
      if( cmd.hasOption("stopForDebug"))
      {  BufferedReader localRdr = new BufferedReader(new InputStreamReader( System.in ));
         System.err.print("Please press return:");

         localRdr.readLine();
      }
   }

   @SuppressWarnings("all")
   public static void printOEObjects(String msg)
   {  if( ! PRINTOEObjects ) return;

      System.err.println("\n" + msg);

      //OENativePtr.PrintLiveObjects();
      System.gc(); System.gc(); System.gc(); System.gc(); System.gc(); System.gc();
      //OENativePtr.PrintLiveObjects();
   }


   private void printRules(PrintStream out) {
      out.println("<html><head><base href='http://research/'/><body>");
      out.println("<h1>Genentech Structure Normalization Rules</h1>");
      out.println("To regenerate this documentation run: 'ant javaDoc'.");
      out.println("<table border='1'><tr><th>Name</th><th>Description</th></tr>\n");

      // print rules that are not implemented as checker because they work on the
      // molfile string
      if( checkForAtomLabelDesc != null )
         out.printf("<tr><td>%s</td><td>%s</td></tr>\n",
               "atomLabelCheck", checkForAtomLabelDesc );

      if( checkForThickBondDesc != null )
         out.printf("<tr><td>%s</td><td>%s</td></tr>\n",
               "thickBondCheck", checkForThickBondDesc );

//      if( checkForWigglyBondDesc != null )
//         out.printf("<tr><td>%s</td><td>%s</td></tr>\n",
//               "wigglyBondCheck", checkForWigglyBondDesc );
//

      // print all other rules
      for( StructureCheckInterface r : rules )  {
         out.printf("<tr><td>%s</td><td>%s</td></tr>\n",
            r.getCheckName(),
            r.getDescriptionHTML());
      }
      out.println("</table></body></html>");
   }

   /**
    *
    * @return may return null if no SFlagChecker was defined.
    */
   public StructureFlag getStructureFlag() {
      if(sflagChecker == null) return null;
      return sflagChecker.getStructureFlag();
   }

   /**
    *
    * @return returns true if the last check found an  input-stereoflag which is
    *          inconsistent with structure.
    */
   public boolean hasStructureFlagError() {
      if(sflagChecker == null) return false;
      return sflagChecker.hasStructureFlagError();
   }

   /**
    * Command line interface to {@link OEStruchk}.
    */
   public static void main2(String [] args) throws ParseException {
      long start = System.currentTimeMillis();
      int nMessages = 0;
      int nErrors = 0;
      int nStruct = 0;

      // create command line Options object
      Options options = new Options();
      Option opt = new Option("f", true, "specify the configuration file name");
      opt.setRequired(true);
      options.addOption(opt);

      CommandLineParser parser = new PosixParser();
      CommandLine cmd = parser.parse( options, args);
      args = cmd.getArgs();

      if(args.length < 1) {
         HelpFormatter formatter = new HelpFormatter();
         formatter.printHelp( "oeStruck", options );
         throw new Error( "missing input file\n");
      }

      String confFile = cmd.getOptionValue("f");

      try {
         oemolistream ifs = new oemolistream(args[0]);
         oemolostream ofs = new oemolostream(args[1]);

         // create OEStruchk from config file
         OEStruchk strchk = new OEStruchk(new File(confFile).toURI().toURL(),
                                          CHECKConfig.TESTAssignStuctFlag, false);

         OEGraphMol mol = new OEGraphMol();
         while ( oechem.OEReadMolecule(ifs , mol ) ) {
            if(! strchk.applyRules(mol, null))
               nErrors++;

            OEGraphMol outMol = strchk.getTransformedMol(null);

            List<Message> msgs = strchk.getStructureMessages(null);
            if(msgs.size() > 0) {
               nMessages += msgs.size();
               StringBuilder sb = new StringBuilder();

               for( Message msg : msgs )
                  sb.append(String.format("%s: %s\n", msg.getLevel(), msg.getText()));

               oechem.OESetSDData(outMol,"errors_oe2", sb.toString());

            }

            oechem.OESetSDData(outMol,"ISM", strchk.getTransformedIsoSmiles(null));
            oechem.OESetSDData(outMol,"SMI", strchk.getTransformedSmiles(null));

            oechem.OEWriteMolecule(ofs, outMol);

            msgs.clear();
         }


         ofs.close();
         ofs.delete();
         ifs.close();
         ifs.delete();
      } catch (Exception e) {
         throw new Error(e);
      }finally {
         System.err.printf("Checked %d structures %d errors, %d messages in %dsec\n",
               nStruct, nErrors, nMessages, (System.currentTimeMillis()-start)/1000);
      }
   }

   public static URL getResourceURL(Class<?> cls, String name)
   {  String path = cls.getName().replace('.','/');
      path = path.subSequence(0,path.lastIndexOf('/')+1) + name;
      return cls.getClassLoader().getResource(path);
   }


   private String countStereoDBondStr() {
      int nStereoDBond = sflagChecker.getNStereoDBond();
      int nStereoDBondSpecified = sflagChecker.getNStereoDBondSpecified();

      return String.format("%s.%s",
         nStereoDBond          == 0 ? "" : Integer.toString(nStereoDBond),
         nStereoDBondSpecified == 0 ? "" : Integer.toString(nStereoDBondSpecified));
   }


   private String countChiralCentersStr() {
      int nChiral          = sflagChecker.getNChiral();
      int nChiralSpecified = sflagChecker.getNChiralSpecified();
      int nNonChiralSp3    = sflagChecker.getNNonChiralStereo();
      int nNonChiralSp3Specified = sflagChecker.getNNonChiralStereoSpecified();

      return String.format("%s.%s.%s.%s",
            nChiral               == 0 ? "" : Integer.toString(nChiral),
            nChiralSpecified      == 0 ? "" : Integer.toString(nChiralSpecified),
            nNonChiralSp3         == 0 ? "" : Integer.toString(nNonChiralSp3),
            nNonChiralSp3Specified== 0 ? "" : Integer.toString(nNonChiralSp3Specified));
   }


   /** this invalidates the {@link OEStruchk} instance */
   public void delete() {
      canonizer.delete();
      oeMol.delete();
      tmpMol.delete();
      structMessages.clear();

      for(StructureCheckInterface rule : rules)
         rule.delete();

      hAtomFunct.delete();
      atomHadExplicitHFunct.delete();
   }

   /** true if current mol is empty and was replaced by "*" */
   public boolean isEmptyMol() {
      return isEmptyMol;
   }


   /** to be called by Checks to make sure the molecule has be properly initilized
    * for all checkes
    */
   public static void initMolForChecks(OEGraphMol mol)
   {  // tag atoms so that we can distinguish them afterwards
      int aNum = 0;
      OEAtomBaseIter aIt = mol.GetAtoms();
      while( aIt.hasNext() )
         aIt.next().SetIntData(ATOMNumTag, aNum++);
      aIt.delete();
   }
}



class AtomBond
{  OEAtomBase otherAtom;
   OEBondBase bd;
   private int otherAtRingCount = -1;

   public AtomBond(OEBondBase bd, OEAtomBase otherAt) {
      this.bd = bd;
      this.otherAtom = otherAt;
   }

   public boolean isAliphaticBridgeHead() {
      if(! bd.IsInRing() ) return false;
      if( otherAtom.GetAtomicNum()!= 6) return false;
      if( otherAtom.GetDegree()   != 4) return false;
      if( otherAtomRings() != 3) return false;

      return true;
   }

   boolean isRingBond() { return bd.IsInRing(); }

   int otherAtomRings () {
      if( otherAtRingCount > -1 ) return otherAtRingCount;

      otherAtRingCount = 0;
      OEBondBaseIter bIt = otherAtom.GetBonds();
      while(bIt.hasNext())
         if(bIt.next().IsInRing()) otherAtRingCount++;
      bIt.delete();

      return otherAtRingCount;
   }
}
