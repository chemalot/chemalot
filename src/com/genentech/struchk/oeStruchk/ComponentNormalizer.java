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

import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import openeye.oechem.OEAtomBase;
import openeye.oechem.OEAtomBaseIter;
import openeye.oechem.OEGraphMol;
import openeye.oechem.OEMolBase;
import openeye.oechem.OEPartPredAtom;
import openeye.oechem.oechem;

import org.jdom.Element;

import com.aestel.io.dataAccess.ConnectionFactory;
import com.aestel.io.dataAccess.ConnectionWrapper;
import com.aestel.utility.Message;
import com.aestel.utility.MessageList;
import com.genentech.oechem.tools.OETools;
import com.genentech.struchk.oeStruchk.OEStruchk.StructureFlag;

/**
 * Structure checking rules that checks for invalid atomic valences.
 *
 * Atoms which have not been specified may have any valence.
 *
 * @author albertgo
 *
 */
public class ComponentNormalizer extends AbstractStructureCheck
              implements ComponentNormalizerInterface, StructureKeeperInterface {

   /** Map from ismi to Salt Object */
   private final Map<String,Salt> saltMap;
   private final Set<String> solvents;
   private final boolean allowMixtures;

   /** Salt code for parent compounds (with no counter ion) */
   private final Salt parentSalt;

   /** stores smiles of salt of last compound */
   private String saltISmi = null;
   /** stores the number of slat molecules of last compound */
   private int saltCount = 0;

   private OEGraphMol keeperSubstance;
   private String substKeeperName;
   private MessageList keeperMsgs;

   ComponentNormalizer(Element tElem) {
      super(tElem);

      allowMixtures = "y".equalsIgnoreCase(tElem.getAttributeValue("allowMixtures"));
      substKeeperName = tElem.getAttributeValue("substanceKeeperName");

      boolean readSaltsFromXML = "xml".equals(tElem.getAttributeValue("saltDefinition"));

      // read salts
      if(readSaltsFromXML)
         saltMap = readXMLSaltMap(tElem);
      else
         saltMap = readDBSaltMap(tElem);

      // parent salt code has empty smiles
      if(saltMap.get("") == null)
         throw new Error("No parent salt code defined!");
      parentSalt = saltMap.get("");

      // read solvents
      solvents =readXMLSolvents(tElem, saltMap);

      if(! checkExample())
         throw new Error( String.format("Example %s did not contain solvents or salts.",
                                 getExampleInput()));
   }

   private static Map<String, Salt> readDBSaltMap(Element elem) {
      Map<String, Salt> sltMap = new HashMap<String, Salt>();
      String sql = elem.getChildText("saltSql");
      ConnectionWrapper con = null;
      Statement stmt = null;
      ResultSet rs = null;
      try {
         con = ConnectionFactory.getDefaultConnection();
         stmt = con.createStatement();
         rs = stmt.executeQuery(sql);
         while(rs.next()) {
            String code = rs.getString(1);
            String smi  = rs.getString(2);
            String name = rs.getString(3);
            if(smi==null) smi ="";

            Salt salt = new Salt(smi, code, name);
            if(sltMap.containsKey(salt.canISmiles))
               throw new Error("Duplicate salts with same smiles in database: "
                              + salt.canISmiles );
            sltMap.put(salt.canISmiles, salt);
         }
      } catch (SQLException e) {
         throw new Error(e);
      }finally {
         try {
            if(rs!=null) rs.close();
            if(stmt!=null) stmt.close();
            if(con!=null)  con.close();
         } catch (SQLException e) {
            throw new Error(e);
         }
      }
      return sltMap;
   }

   private static Map<String, Salt> readXMLSaltMap(Element tElem) {
      Map<String, Salt> sltMap = new HashMap<String, Salt>();
      for(Object saltElementO :  tElem.getChildren("salt")) {
         Element saltElement = (Element)saltElementO;
         String smi = saltElement.getAttributeValue("smiles");
         String code   = saltElement.getAttributeValue("code");

         Salt salt = new Salt(smi, code, code);
         sltMap.put(salt.canISmiles, salt);
      }
      return sltMap;
   }

   private static Set<String> readXMLSolvents(Element tElem, Map<String, Salt> saltMap) {
      Set<String>solvents = new HashSet<String>();
      OEGraphMol mol = new OEGraphMol();
      for(Object solvElementO :  tElem.getChildren("solvent")) {
         Element solvElement = (Element)solvElementO;
         String smi = solvElement.getAttributeValue("smiles");

         OETools.smiToMol(mol, smi);
         if(mol.NumAtoms() == 0)
            throw new Error(String.format("Invalid solvent smiles: %s.", smi));

         String iSmi = OETools.molToCanSmi(mol, true);
//         This can only be activated when the solvent salcodes are removed from table after migration in 2009
//         if( saltMap.containsKey(iSmi) )
//         {  if( !Settings.ISProduction )
//               throw new Error("Trying to load solvent for a salt: " + iSmi);
//            else
//               System.err.println("\n\nTrying to load solvent for a salt: " + iSmi + "\n\n");
//         }else
//         {
               solvents.add(iSmi);
//         }

         mol.Clear();
      }
      mol.delete();

      return solvents;
   }

   @Override
   public boolean checkStructure(OEGraphMol in, StructureFlag inStereo,
                                 MessageList msgs) {
      saltISmi = null;
      saltCount = 0;

      int[] parts = new int[in.GetMaxAtomIdx()];
      int nComponents = oechem.OEDetermineComponents(in, parts);

      if(nComponents == 1) {
         if(substKeeperName != null) keepSubstance(in, msgs);
         return true;
      }

      MolComponent[] components = new MolComponent[nComponents];
      OEPartPredAtom pred = new OEPartPredAtom(parts);

      // get components from in
      for (int i = 1; i <= nComponents; i++) {
        pred.SelectPart(i);
        OEGraphMol partMol = new OEGraphMol();
        oechem.OESubsetMol(partMol, in, pred);

        components[i-1] = new MolComponent(i, partMol);
      }
      pred.delete();

      // sort by smiles length so that we can start with the shortest to remove
      Arrays.sort(components);

      MolComponent last = null;
      // remove identical components
      for(int i=0; i<components.length; i++) {
         MolComponent mc = components[i];
         if(mc == null) continue;

         if(last != null && mc.getCanISmi().equals(last.getCanISmi())) {
            removeComponent(in, parts, mc.getId());
            last.incrementOccurenceCount();

            mc.delete();
            components[i] = null;
            nComponents--;
            continue;
         }
         last = mc;
      }

      // remove solvents
      for(int i=0; i<components.length; i++) {
         if(components[i] == null) continue;

         // remove solvent
         if(solvents.contains(components[i].getCanISmi()) && nComponents > 1) {
            if(components[i] == null) continue;

            removeComponent(in, parts, components[i].getId());
            msgs.addMessage(new Message("Removed solvent: " + components[i].getCanISmi(),
                  Message.Level.COMMENT, null));
            components[i].delete();
            components[i] = null;
            nComponents--;
            continue;
         }
      }

      // keep substance molecule with single copy of counter ion
      if(substKeeperName != null) keepSubstance(in, msgs);

      // make "*" component last so that other salts get replaced first
      // "*" can be either the smiles of a nullstructure or the salt signifying a micture of salts
      for(int i=0; i<components.length; i++) {
         if( components[i] != null && "*".equals(components[i].getCanISmi()) )
         {  MolComponent dummy = components[components.length-1];
            components[components.length-1] = components[i];
            components[i] = dummy;
            break;
         }
      }

      // remove salts
      for(int i=0; i<components.length; i++) {
         if(components[i] == null) continue;

         if(saltMap.containsKey(components[i].getCanISmi()) && nComponents > 1 ) {
            // remove salt
            removeComponent(in, parts, components[i].getId());
            msgs.addMessage(new Message("Removed counter ion: " + components[i].getCanISmi(),
                  Message.Level.COMMENT, null));

            if(saltISmi != null)
            {  msgs.addMessage(new Message(
                  String.format("Mixed counter ion found: %s and %s",
                  saltISmi, components[i].getCanISmi()), Message.Level.WARNING, null));
               saltISmi = "*";
               saltCount= 1;
            } else
            {  saltISmi = components[i].getCanISmi();
               saltCount = components[i].getOccurenceCount();
            }

            components[i].delete();
            components[i] = null;
            nComponents--;
         }
      }

      if("*".equals(saltISmi) && substKeeperName != null )
      {  // Mixed salt, replace counter ions with *
         replaceKeepSubstanceMixedSalt(in, msgs);
      }

      if(allowMixtures || inStereo == StructureFlag.UNCERTAINStructure ) {
         deleteComponents(components);    // clean memory
         return true;
      }

      // check that left components have same molecular formula (are isomers)
      String lastMF = null;
      if(nComponents > 1) {
         for(int i=0; i<components.length; i++) {
            if(components[i] == null) continue;

            String mf = oechem.OEMolecularFormula(components[i].getMol());
            if(! mf.equals(lastMF) && lastMF != null) {
               msgs.addMessage(new Message("Structure has multiple components.", Message.Level.ERROR, null));
               deleteComponents(components);
               return false;
            }
            components[i].delete();
            lastMF = mf;
         }
      }

      deleteComponents(components);
      return true;
   }

   private void deleteComponents(MolComponent[] components) {
      for(int i=0; i<components.length; i++) {
         if(components[i] == null) continue;

         components[i].delete();
      }
   }

   @Override
   public String getSaltCode()
   {  if(saltISmi == null) return getParentSaltCode();

      return saltMap.get(saltISmi).code;
   }

   @Override
   public String getSaltName()
   {  if(saltISmi == null) return getParentSaltName();

      return saltMap.get(saltISmi).name;
   }

   /** return the number of repetitions of the salt molecule in the last Structure */
   @Override
   public int getSaltCount()
   {  return saltCount;
   }

   @Override
   public String getSaltMF() {
      if(saltISmi == null) return "";

      return saltMap.get(saltISmi).mf;
   }

   @Override
   public String getSaltMW() {
      if(saltISmi == null) return "0";

      return saltMap.get(saltISmi).mw;
   }


   /**
    * remove atoms from mol which are part of component id.
    * @param parts a map from atom.getIdx() to component id
    */
   private void removeComponent(OEGraphMol in, int[] parts, int id) {
      OEAtomBaseIter aIt = in.GetAtoms();
      while( aIt.hasNext() ) {
         OEAtomBase at = aIt.next();
         if(parts[at.GetIdx()] == id) in.DeleteAtom(at);
      }
      aIt.delete();
   }

   @Override
   public String getParentSaltCode() {
      return parentSalt.code;
   }

   public String getParentSaltName() {
      return parentSalt.name;
   }

   @Override
   public void reset() {
      saltISmi = null;
   }

   @Override
   public void delete() {
      if( keeperSubstance != null) keeperSubstance.delete();
      keeperSubstance = null;
   }


   ////////////////////////////// {@link StructureKeeperInterface}

   private void keepSubstance(OEGraphMol in, MessageList msgs) {
      keeperSubstance = new OEGraphMol(in); // keep a copy so nobody can change it

      oechem.OESuppressHydrogens(keeperSubstance,false,false,true);
      oechem.OEPerceiveChiral(keeperSubstance);

      this.keeperMsgs = new MessageList(msgs);
   }


   private static final OEMolBase WILDCARDMolecule = OETools.smiToMol(new OEGraphMol(), "*");

   private void replaceKeepSubstanceMixedSalt(OEGraphMol in, MessageList msgs) {
      keeperSubstance.Clear();
      oechem.OEAddMols(keeperSubstance, in);

      OETools.combineStructures(keeperSubstance, WILDCARDMolecule);

      oechem.OESuppressHydrogens(keeperSubstance,false,false,true);
      oechem.OEPerceiveChiral(keeperSubstance);

      this.keeperMsgs = new MessageList(msgs);
   }

   @Override
   public MessageList getStructureMessages() {
      return keeperMsgs;
   }

   @Override
   public String getKeeperName() {
      return substKeeperName;
   }

   @Override
   public OEGraphMol getMolecule() {
      return keeperSubstance;
   }
}
