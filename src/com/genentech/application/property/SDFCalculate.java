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
package com.genentech.application.property;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.URL;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Vector;
import openeye.oechem.*;
import org.apache.commons.cli.*;
import org.jdom.Element;

import com.aestel.io.XMLUtil;
import com.aestel.math.StepFunction;

/*main class to start calculating TPSA and other 2D properties */
public class SDFCalculate {

   private static final String CLOGD74_Tag = "cgLogD7.4";

   private final boolean countP;
   private final boolean countS;
   private final String cLogPTag;
   private final Map<String, String> smartsMap;
   private final Outputter out;
   private static final String protonateCNOS="[O,S,#7,#6;-1;!$(*[*+]):1]>>[*+0:1][H]";
   private static final String deprotonateNP="[#7,#15;+1:1][H]>>[*;+0:1]";

   private static  Map<String,String> propsMap = createMap();

   private static Map<String, String> createMap() {
      Map<String, String> result = new LinkedHashMap<String, String>();
      result.put("all","all properties");
      result.put("CNS_MPO","CNS_MPO score, requires cLogP, pKa_MB, cgLogD7.4 ");
      result.put("Charge","ionized by MoKa at pH7.4");
      result.put("cIC50atLE0.3","calculated IC50 at LE 0.3");
      result.put("cIC50atLE0.35","calculated IC50 at LE 0.35");
      result.put("cIC50atLE0.4","calculated IC50 at LE 0.4");
      result.put("Heavy_Atoms","Num of heavy atoms");
      result.put("H_polar","ionized by Moka at pH7.4");
      result.put("MW","neutralized by OpenEye's OEchem");
      result.put("N+O","number of N + number of O");
      result.put("NH+OH","number of NH + OH, NH2 counts as 2");
      result.put("Rings","number of ring systems");
      result.put("RotBonds","number of ratatable bonds");
      result.put("RO5"," number of Lipinski's rule of five violations");
      result.put("TPSA","topological polar surface area");
      result.put("AromaticFraction","number of aromatic atoms / number of heavy atoms");
      result.put("CarboAromaticFraction","number of aromatic carbon atoms / number of heavy atoms");
      result.put("AromaticRings","number of 5, 6, or 7 member aromatic rings");
      result.put("CarboAromaticRings","number of 5, 6, or 7 member carboaromatic rings");
      result.put("HeteroAromaticRings","number of 5, 6, or 7 member heteroaromatic rings");
      result.put("AliphaticRings","number of 3,4,5,6,7,or 8 member aliphatic rings ");
      result.put("CarboAliphaticRings","number of 3,4,5,6,7,or 8 member carboaliphatic rings ");
      result.put("HeteroAliphaticRings","number of 3,4,5,6,7,or 8 member heteroaliphatic rings ");
      result.put("Solubility_Index","cgLogD7.4 + AromaticRings");
      result.put("Csp3","number of sp3 carbons");
      result.put("Csp3Fraction","Fraction of sp3 carbons");
      result.put("NonSp3Fraction","fraction of non-sp3 atoms");
      result.put("LargestRingSize","atom count of largest ring");      
      result.put("TotalAtoms","atom count inlcuding any explicit hydrogens");      
      return Collections.unmodifiableMap(result);
  }

   /**
    * @param countP if true phosphorus is considered for tpsa
    * @param countS if true sulfur is considered for tpsa
    * @param cLogPTag name of tag in input containing logP value for rule of five violations
    * @param smartsFile xml file with smarts definition for atom types
    * @throws IOException
    */
   public SDFCalculate(String outFName, String cLogPTag, boolean countP, boolean countS,
            Element root) throws IOException
   {  this.cLogPTag = cLogPTag == null ? "" : cLogPTag;
      this.countS = countS;
      this.countP = countP;
      parseXML myXML = new parseXML();
      this.smartsMap =myXML.parse(root);
      out = new OEOutputter(outFName);
         }

   private void calcProperties(OEGraphMol mol,Vector<String> propsList) {
      TPSA myTPSA = new TPSA();
      int numHeavy=0;
      int hPolar_neu =0;
      double tpsa=0.0;
      double  mw =0;
      int LipinskiHBA=0, LipinskiHBD=0;
      int aromaticRingCount=0, carboAromaticRingCount=0, aliphaticRingCount=0, carboAliphaticRingCount=0;
      final double R=8.314510, J2kcal=0.001*(1./4.184), T=300.0, RT=R*J2kcal*T;

      // assign OEChem default aromatic model and hybridazation
      oechem.OEAssignAromaticFlags(mol);
      oechem.OEAssignHybridization(mol);
      // transform molecule to protonateCNOS and deprotonateNP
      // explicit hydrogens are added to neutralMol
      OEGraphMol neutralMol = neutralizeMol(mol);
      // convert neutralMol to smiles, print smiles
      //String cansmi = oechem.OECreateSmiString(neutralMol);
      //System.err.println("SMILES: " + cansmi);

      //TPSA RotBonds H_polar N O NH OH Heavy_Atoms Rings MW RO5 cLogP_read Charge
      if (propsList.contains("TPSA") || propsList.contains("CNS_MPO") || propsList.contains("all")) {
         tpsa = myTPSA.calculateTPSA(neutralMol, countP, countS);
         if (tpsa <=0.0) tpsa = 0.0;
         oechem.OESetSDData(mol, "TPSA", Integer.toString((int)tpsa));
      }
      if (propsList.contains("Charge") || propsList.contains("all") ) {
         int netCharge = oechem.OENetCharge(mol);
         oechem.OESetSDData(mol, "Charge",     Integer.toString(netCharge));
      }
      if (propsList.contains("RotBonds") ||propsList.contains("all") ) {
         /* use SMARTS to calculate the following properties */
         int singleBonds = SmartsSearch.search(mol, smartsMap.get("SingleBond"));
         int amideCount = SmartsSearch.search(mol, smartsMap.get("Amide"));
         int tripleBondCount = SmartsSearch.search(mol,smartsMap.get("TripleBondAtom"));
         int rotBonds = singleBonds - amideCount - tripleBondCount;
         oechem.OESetSDData(mol, "RotBonds",   Integer.toString(rotBonds));
      }
      if (propsList.contains("Rings") || propsList.contains("all") ) {
         //oechem.OEFindRingAtomsAndBonds(mol); // doesn't seem to be necessary
         int numRings = oechem.OEDetermineRingSystems(mol, new int[mol.GetMaxAtomIdx()]);
         oechem.OESetSDData(mol, "Rings",      Integer.toString(numRings));
      }
      if (propsList.contains("LargestRingSize") || propsList.contains("all") ) {         
         int maxRingSize     = 0;          
         int currentRingSize = 0;
         for (OEAtomBase atom : mol.GetAtoms()) {
            currentRingSize = oechem.OEAtomGetSmallestRingSize(atom);
            if (currentRingSize > maxRingSize) {
               maxRingSize = currentRingSize;
            }
         }
         oechem.OESetSDData(mol, "LargestRingSize",      Integer.toString(maxRingSize));               
      }
      if (propsList.contains("Heavy_Atoms") || propsList.contains("cIC50atLE0.3") || propsList.contains("all") ) {
         numHeavy = SmartsSearch.search(mol, smartsMap.get("HeavyAtom"));
         oechem.OESetSDData(mol, "Heavy_Atoms",Integer.toString(numHeavy));
      }
      if (propsList.contains("cIC50atLE0.3") || propsList.contains("all") ) {
         double cIC50atLE0_3=Math.exp((-1.0*0.3*numHeavy)/RT)*1000000; // multiply by 1000000 to convert to uM units
         oechem.OESetSDData(mol, "cIC50atLE0.3",String.format("%.5f", cIC50atLE0_3));
      }
      if (propsList.contains("cIC50atLE0.35") || propsList.contains("all") ) {
         double cIC50atLE0_3=Math.exp((-1.0*0.35*numHeavy)/RT)*1000000; // multiply by 1000000 to convert to uM units
         oechem.OESetSDData(mol, "cIC50atLE0.35",String.format("%.5f", cIC50atLE0_3));
      }
      if (propsList.contains("cIC50atLE0.4") || propsList.contains("all") ) {
         double cIC50atLE0_3=Math.exp((-1.0*0.4*numHeavy)/RT)*1000000; // multiply by 1000000 to convert to uM units
         oechem.OESetSDData(mol, "cIC50atLE0.4",String.format("%.5f", cIC50atLE0_3));
      }
      if (propsList.contains("NH+OH") || propsList.contains("RO5") || propsList.contains("CNS_MPO") || propsList.contains("all") ) {
         hPolar_neu = countHPolar(neutralMol);
         LipinskiHBD = hPolar_neu;
         oechem.OESetSDData(mol, "NH+OH",      Integer.toString(hPolar_neu));
      }
      if (propsList.contains("H_polar") || propsList.contains("all") ) {
         int hPolar = countHPolar(mol);
         oechem.OESetSDData(mol, "H_polar",    Integer.toString(hPolar));
      }
      if (propsList.contains("N+O") || propsList.contains("RO5") || propsList.contains("all") ) {
         int nCount = SmartsSearch.search(mol, smartsMap.get("NCount"));
         int oCount = SmartsSearch.search(mol, smartsMap.get("OCount"));
         LipinskiHBA = nCount + oCount;
         oechem.OESetSDData(mol, "N+O",        String.format("%d", LipinskiHBA));
      }
      if (propsList.contains("MW") || propsList.contains("CNS_MPO") || propsList.contains("RO5") || propsList.contains("all") ) {
         mw =   oechem.OECalculateMolecularWeight(neutralMol, true);
         oechem.OESetSDData(mol, "MW",         String.format("%.2f ", mw));
      }
      if (propsList.contains("AromaticRings") || propsList.contains("HeteroAromaticRings")
       || propsList.contains("Solubility_Index")  || propsList.contains("all") ) {
         int aCount = SmartsSearch.search(mol, smartsMap.get("Aromatic5Rings"));
         aCount    += SmartsSearch.search(mol, smartsMap.get("Aromatic6Rings"));
         aCount    += SmartsSearch.search(mol, smartsMap.get("Aromatic7Rings"));
         aromaticRingCount = aCount;
         oechem.OESetSDData(mol, "AromaticRings",         String.format("%d", aCount ));
      }
      if (propsList.contains("CarboAromaticRings") || propsList.contains("HeteroAromaticRings") || propsList.contains("all") ) {
         int caCount = SmartsSearch.search(mol, smartsMap.get("CarboAromatic5Rings"));
         caCount    += SmartsSearch.search(mol, smartsMap.get("CarboAromatic6Rings"));
         caCount    += SmartsSearch.search(mol, smartsMap.get("CarboAromatic7Rings"));
         carboAromaticRingCount = caCount;
         oechem.OESetSDData(mol, "CarboAromaticRings",     String.format("%d", caCount ));
      }
      if (propsList.contains("HeteroAromaticRings") || propsList.contains("all") ) {
         int haCount = aromaticRingCount - carboAromaticRingCount;
         oechem.OESetSDData(mol, "HeteroAromaticRings",     String.format("%d", haCount ));
      }
      if (propsList.contains("AliphaticRings") || propsList.contains("HeteroAliphaticRings") || propsList.contains("all") ) {
         int aCount = SmartsSearch.searchAliphaticRings(mol, smartsMap.get("Aliphatic3Rings"));
         aCount    += SmartsSearch.searchAliphaticRings(mol, smartsMap.get("Aliphatic4Rings"));
         aCount    += SmartsSearch.searchAliphaticRings(mol, smartsMap.get("Aliphatic5Rings"));
         aCount    += SmartsSearch.searchAliphaticRings(mol, smartsMap.get("Aliphatic6Rings"));
         aCount    += SmartsSearch.searchAliphaticRings(mol, smartsMap.get("Aliphatic7Rings"));
         aCount    += SmartsSearch.searchAliphaticRings(mol, smartsMap.get("Aliphatic8Rings"));
         aliphaticRingCount = aCount;
         oechem.OESetSDData(mol, "AliphaticRings",         String.format("%d", aCount ));
      }

      if (propsList.contains("CarboAliphaticRings") || propsList.contains("HeteroAliphaticRings") || propsList.contains("all") ) {
         int aCount = SmartsSearch.searchAliphaticRings(mol, smartsMap.get("CarboAliphatic3Rings"));
         aCount    += SmartsSearch.searchAliphaticRings(mol, smartsMap.get("CarboAliphatic4Rings"));
         aCount    += SmartsSearch.searchAliphaticRings(mol, smartsMap.get("CarboAliphatic5Rings"));
         aCount    += SmartsSearch.searchAliphaticRings(mol, smartsMap.get("CarboAliphatic6Rings"));
         aCount    += SmartsSearch.searchAliphaticRings(mol, smartsMap.get("CarboAliphatic7Rings"));
         aCount    += SmartsSearch.searchAliphaticRings(mol, smartsMap.get("CarboAliphatic8Rings"));
         carboAliphaticRingCount = aCount;
         oechem.OESetSDData(mol, "CarboAliphaticRings",         String.format("%d", aCount ));
      }
      if (propsList.contains("HeteroAliphaticRings") || propsList.contains("all") ) {
         int count = aliphaticRingCount - carboAliphaticRingCount;
         oechem.OESetSDData(mol, "HeteroAliphaticRings",         String.format("%d", count ));
      }

      if (propsList.contains("AromaticFraction") || propsList.contains("all") ) {
         int aCount = SmartsSearch.search(mol, smartsMap.get("AromaticAtom"));
         int heavyCount = SmartsSearch.search(mol, smartsMap.get("HeavyAtom"));
         double aromaticFraction = aCount / (heavyCount + 0.001); // so you don't divide by zero
         oechem.OESetSDData(mol, "AromaticFraction",         String.format("%.2f", aromaticFraction ));
      }
      if (propsList.contains("CarboAromaticFraction") || propsList.contains("all") ) {
         int cCount = SmartsSearch.search(mol, smartsMap.get("AromaticCarbonAtom"));
         int heavyCount = SmartsSearch.search(mol, smartsMap.get("HeavyAtom"));
         double carboAromaticFraction = cCount / (heavyCount + 0.001); // so you don't divide by zero
         oechem.OESetSDData(mol, "CarboAromaticFraction",         String.format("%.2f", carboAromaticFraction ));
      }
      
      if (propsList.contains("NonSp3Fraction") || propsList.contains("all") ) {
         int nonSp3Count = SmartsSearch.search(mol, smartsMap.get("sp3"));
         int heavyCount = SmartsSearch.search(mol, smartsMap.get("HeavyAtom"));
         double nonSp3Fraction = (heavyCount - nonSp3Count) / (heavyCount + 0.001); // so you don't divide by zero
         oechem.OESetSDData(mol, "NonSp3Fraction",         String.format("%.2f", nonSp3Fraction ));
      }
      
      if (propsList.contains("CarboAromaticRings") || propsList.contains("HetereoAromaticRings") || propsList.contains("all") ) {
         int caCount = SmartsSearch.search(mol, smartsMap.get("CarboAromatic5Rings"));
         caCount    += SmartsSearch.search(mol, smartsMap.get("CarboAromatic6Rings"));
         caCount    += SmartsSearch.search(mol, smartsMap.get("CarboAromatic7Rings"));
         carboAromaticRingCount = caCount;
         oechem.OESetSDData(mol, "CarboAromaticRings",     String.format("%d", caCount ));
      }
      
      if (propsList.contains("Csp3") || propsList.contains("all") ) {
         int count = SmartsSearch.search(mol, smartsMap.get("Csp3"));
         oechem.OESetSDData(mol, "Csp3", String.format("%d", count ));

         count = SmartsSearch.search(mol, smartsMap.get("CSsp3"));
         oechem.OESetSDData(mol, "CSsp3", String.format("%d", count ));

         count = SmartsSearch.search(mol, smartsMap.get("CS2sp3"));
         oechem.OESetSDData(mol, "CS2sp3", String.format("%d", count ));
      }      

      if (propsList.contains("Csp3Fraction") || propsList.contains("all") ) {
         int heavyCount = SmartsSearch.search(mol, smartsMap.get("HeavyAtom"));
         int count = SmartsSearch.search(mol, smartsMap.get("Csp3"));
         double Csp3Fraction = (count / (heavyCount + 0.001)); // so you don't divide by zero
         oechem.OESetSDData(mol, "Csp3Fraction", String.format("%.2f", Csp3Fraction ));
      }
      if (propsList.contains("TotalAtoms") || propsList.contains("all") ) {
         OEAtomBaseIter atIt = mol.GetAtoms();
         int implH = 0;
         while(atIt.hasNext())
            implH += atIt.next().GetImplicitHCount();
         atIt.delete();
         oechem.OESetSDData(mol, "TotalAtoms", String.format("%d", mol.NumAtoms()+implH ));
      }

      
      // get cLogP for RO5 and CNS_MPO calculation
      double cLogP =Double.NaN;
      if (propsList.contains("RO5") || propsList.contains("CNS_MPO") || propsList.contains("all") ) {
         if (cLogPTag != null ){
            if (oechem.OEHasSDData(mol, cLogPTag)) {
               try {
                  cLogP = Double.valueOf(oechem.OEGetSDData(mol, cLogPTag));
               }catch (Exception e){
                  System.err.println(mol.GetTitle() + " " + oechem.OEGetSDData(mol, cLogPTag));
               }
            }
         }
      }

      //get cLogD74 for CNS_MPO and Solubility_index calculation
      double cLogD74=Double.NaN;
      if (propsList.contains("CNS_MPO") || propsList.contains("Solubility_Index") || propsList.contains("all") ) {
         if (oechem.OEHasSDData(mol, CLOGD74_Tag)) {
            try {
               cLogD74 = Double.valueOf(oechem.OEGetSDData(mol, CLOGD74_Tag));
            }catch (Exception e){
               System.err.println("mol: " + mol.GetTitle() + " " + CLOGD74_Tag+ ": " + oechem.OEGetSDData(mol, CLOGD74_Tag));
               cLogD74 = Double.NaN;
            }
         }
      }

      if (propsList.contains("Solubility_Index") || propsList.contains("all") ) {
         if ( Double.isNaN(cLogD74) ) { // do Solubility_Index if cLogD74 is not found
            System.err.println("mol: "+mol.GetTitle() + " Solubility_Index error, cLogD74 is not defined" );
         }else {
            double solIdx = cLogD74 + aromaticRingCount;
            oechem.OESetSDData(mol, "Solubility_Index",     String.format("%.1f", solIdx ));
         }
      }

      if (propsList.contains("RO5") || propsList.contains("all") ) {
         /*rule of Five violations, include cLogP if it is provided */
         int ruleOf5Violation =0;
         if (mw > 500.0) ruleOf5Violation++;
         if (LipinskiHBA > 10 )ruleOf5Violation++;
         if (LipinskiHBD > 5 ) ruleOf5Violation++;
         if (cLogP > 5.0) ruleOf5Violation++;
         if ( Double.isNaN(cLogP) ) { // do not calculate RO5 violation if cLogP is not found
            System.err.println("mol: "+mol.GetTitle() + " RO5 error: cLogP is not defined" );
         }else {
            oechem.OESetSDData(mol, "RO5",        Integer.toString(ruleOf5Violation));
         }
      }
      if (propsList.contains("CNS_MPO") || propsList.contains("all") ) {
         double c_pKa_MB=Double.NaN;
         String c_pKa_MB_Tag="c_pKa_MB";
         // if c_pKa_MB does not exist, add one to CNS_MPO score
         if (oechem.OEHasSDData(mol, c_pKa_MB_Tag)) {
            try {
               c_pKa_MB = Double.valueOf(oechem.OEGetSDData(mol, c_pKa_MB_Tag));
            }catch (Exception e){
               System.err.println("mol: "+mol.GetTitle() + " c_pKa_MB: " + oechem.OEGetSDData(mol, c_pKa_MB_Tag));
               c_pKa_MB = Double.NaN;
            }
         }else {
            c_pKa_MB=0.0; // the step function will assign a CNS_MPO_pKa score of 1 for this pka_MB value
         }

         if ( Double.isNaN(cLogP) || Double.isNaN(cLogD74) || Double.isNaN(c_pKa_MB) ){
            //oechem.OESetSDData(mol, "CNS_MPO_score", "error");
            System.err.println("mol: "+mol.GetTitle() + " CNS_MPO_score error" );
         }else {
            double CNS_MPO_cLogP = StepFunction.linearScore(5,3, cLogP);
            double CNS_MPO_cLogD = StepFunction.linearScore(4,2, cLogD74);
            double CNS_MPO_MW    = StepFunction.linearScore(500,360, mw);
            double CNS_MPO_TPSA  = StepFunction.humpScore(20,40,90,120, tpsa);
            double CNS_MPO_HBD   = StepFunction.linearScore(3.5,0.5, hPolar_neu);
            double CNS_MPO_pKa   = StepFunction.linearScore(10,8, c_pKa_MB);

            double CNS_MPO = CNS_MPO_cLogP + CNS_MPO_cLogD + CNS_MPO_MW + CNS_MPO_TPSA + CNS_MPO_HBD + CNS_MPO_pKa;
            oechem.OESetSDData(mol, "CNS_MPO_score",  String.format("%.2f", CNS_MPO));

            // store MPO scores for the individual parts
            oechem.OESetSDData(mol, "CNS_MPO_cLogP",  String.format("%.2f", CNS_MPO_cLogP));
            oechem.OESetSDData(mol, "CNS_MPO_cLogD",  String.format("%.2f", CNS_MPO_cLogD));
            oechem.OESetSDData(mol, "CNS_MPO_MW",  String.format("%.2f", CNS_MPO_MW));
            oechem.OESetSDData(mol, "CNS_MPO_TPSA",  String.format("%.2f", CNS_MPO_TPSA));
            oechem.OESetSDData(mol, "CNS_MPO_HBD",  String.format("%.2f", CNS_MPO_HBD));
            oechem.OESetSDData(mol, "CNS_MPO_pKa",  String.format("%.2f", CNS_MPO_pKa));
         }
      }

   }

   private static OEGraphMol neutralizeMol(OEGraphMol mol) throws Error
   {
      OEGraphMol neutralMol = new OEGraphMol(mol);
      oechem.OEAddExplicitHydrogens(neutralMol);
      OEUniMolecularRxn transform = new OEUniMolecularRxn(protonateCNOS);
      if(! transform.IsValid()) throw new Error("Invalid Smirks " + protonateCNOS);

      transform.constCall(neutralMol);
      transform = new OEUniMolecularRxn(deprotonateNP);
      if(! transform.IsValid()) throw new Error("Invalid Smirks " + deprotonateNP);
      transform.constCall(neutralMol);
      return neutralMol;
   }

   private static int countHPolar(OEGraphMol mol) {
      int HPolar=0;
      OEGraphMol temp = new OEGraphMol (mol);
      oechem.OEAddExplicitHydrogens(temp);

      for(OEAtomBaseIter iter = temp.GetAtoms(new OEIsPolarHydrogen()); iter.hasNext(); ){
         iter.next();
         HPolar++;
      }
      return HPolar;
   }

   /**
    * Calculate properties of all molecules in a file
    *
    * @param inFile
    * @param propsList
    * @param outFile
    * @throws IOException
    */
   public void calcProperties(String inFile, Vector<String> propsList) throws IOException {
      OEGraphMol mol = new OEGraphMol();
      oemolistream ifs;

         ifs = new oemolistream(inFile);
         mol.Clear();
         while (oechem.OEReadMolecule(ifs, mol) ) {
            calcProperties(mol, propsList);
            out.output(mol);
         }
      out.close();
   }

   public static void main(String args[]){
      String usage = "java SDFCalculate [options] <list of space separated properties>\n";

      Options options = new Options();
      // add  options
      options.addOption("TPSA_P", false, "Count phosphorus atoms, default is false. (optional)");
      options.addOption("TPSA_S", false, "Count sulfur atoms, default is false. (optional)");
      options.addOption("cLogP", true, "SDtag where cLogP is stored, default is cLogP (optional)");
      options.addOption("in", true, "inFile in OE formats: Ex: a.sdf or .sdf");
      options.addOption("out", true, "outputfile in OE formats. Ex: a.sdf or .sdf ");

      try {
         boolean countS = false;
         boolean countP = false;

         // append list of valid properties and their descriptions to the usage statement
         Iterator<Entry<String, String>> i = propsMap.entrySet().iterator();
         while(i.hasNext()) {
            Map.Entry<String,String> me = i.next();
            usage = usage + me.getKey() + ":\t" + me.getValue() +"\n";
         }

         CommandLineParser parser = new PosixParser();
         CommandLine cmd = parser.parse(options, args);

         if (cmd.hasOption("TPSA_P"))   countP = true;
         if (cmd.hasOption("TPSA_S")) countS = true;

         // get list of properties
         Vector<String> propsList = new Vector<String>(Arrays.asList(cmd.getArgs()));
         if (propsList.isEmpty()) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp( usage, options );
            System.exit(1);
         }
         //make sure list of requested pros are valid props
         for (String p : propsList){
            if ( ! propsMap.containsKey(p) ){
               System.err.println(p + " is not a valid property.");
               HelpFormatter formatter = new HelpFormatter();
               formatter.printHelp( usage, options );
               System.exit(1);
            }
         }

         // get cLogP SD label tag option value
         String cLogPTag ="cLogP";

         if (cmd.hasOption("cLogP")){
            cLogPTag = cmd.getOptionValue("cLogP");
         }

         String inFile = cmd.getOptionValue("in");
         if (inFile == null) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp( usage, options );
            System.exit(1);
         }

         String outFile = cmd.getOptionValue("out");
         if (outFile == null) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp( usage, options );
            System.exit(1);
         }

         String filename = "smarts.xml";
         URL url = SDFCalculate.class.getResource(filename);
         Element root = XMLUtil.getRootElement(url, false);

         SDFCalculate test = new SDFCalculate(outFile, cLogPTag, countP, countS, root);
         test.calcProperties(inFile, propsList);

      } catch (ParseException e)
      {  HelpFormatter formatter = new HelpFormatter();
         formatter.printHelp( usage, options );
         System.exit(1);
      } catch (Exception e)
      {  e.printStackTrace();
      }
   }
}

interface Outputter
{  void output(OEMolBase mol) throws IOException;
   void close();
}

class OEOutputter implements Outputter
{  oemolostream ofs;

   OEOutputter(String fName) throws FileNotFoundException
   {  ofs = new oemolostream(fName);
   }

   @Override
   public void output(OEMolBase mol)
   {  oechem.OEWriteMolecule(ofs, mol);
   }

   @Override
   public void close()
   {  ofs.close();
   }
}
