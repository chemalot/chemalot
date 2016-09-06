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
import openeye.oequacpac.OETautomerOptions;
import openeye.oequacpac.oequacpac;

import org.jdom.Element;

import com.aestel.utility.Message;
import com.aestel.utility.MessageList;
import com.genentech.oechem.tools.OETools;
import com.genentech.struchk.oeStruchk.OEStruchk.StructureFlag;


/**
 * Canonicalize a structure using quacpac.
 *
 * @author albertgo
 *
 */
public class TautomerStandardizer extends AbstractStructureCheck {
   static final int MAX_TAUTOMER_ENUMERATION = 256;   //maximum number of tautomers to enumerate
   private final TautomerMolEvaluator tEvaluator;
   private final OETautomerOptions tautomerOptions;

   /**
    * Create a Transformer from the xml element.
    */
   TautomerStandardizer(Element tElem) {
      super(tElem);

      tautomerOptions = new OETautomerOptions();
      tautomerOptions.SetMaxTautomersGenerated(MAX_TAUTOMER_ENUMERATION);
      tautomerOptions.SetApplyWarts(false);
      tautomerOptions.SetCarbonHybridization(false);
      tautomerOptions.SetLevel(0);
      tautomerOptions.SetMaxTautomericAtoms(70);
      tautomerOptions.SetMaxZoneSize(35);
      tautomerOptions.SetRankTautomers(false);
      tautomerOptions.SetSaveStereo(true);

      tEvaluator = new TautomerMolEvaluator();

      if(! checkExample())
         throw new Error(
               String.format("Example %s was not transformed by TautomerStandardizer",
                             getExampleInput()));
   }


   @Override
   public boolean checkStructure(OEGraphMol in, StructureFlag inStereo,
                                 MessageList msgs) {
      OEGraphMol tmpMol = new OEGraphMol(in);// work on copy Enumerate changes mol

      String preSmi = OETools.molToCanSmi(tmpMol, true);

      // we could try to set save Stereo to true and rank to true to get better tautomers
      // this would allow us to remove code from TautomerMolEvaluator,
      // however for now lets keep as few changes as possible.

      int nTaut = 0;
      tEvaluator.setupNewMol(tmpMol);
      for (OEMolBase taut : oequacpac.OEEnumerateTautomers(tmpMol, tautomerOptions))
      {  //System.err.printf("%s\t%s\n", oechem.OECreateSmiString(taut,OESMILESFlag.ISOMERIC), OETools.molToCanSmi(taut, true));
         tEvaluator.evaluate(taut);
         taut.delete();
         nTaut++;
      }

      String postSmi = preSmi;
      if(nTaut > 0) {
         OEMolBase tautMol = tEvaluator.GetBest();
         postSmi = OETools.molToCanSmi(tautMol, true);

         if(! preSmi.equals(postSmi)) {
            msgs.addMessage(new Message("Structure converted to default tautomer.",
                  Message.Level.COMMENT , in));
            in.Clear();          // copy transformed into in
            oechem.OEAddMols(in, tautMol);
         }

      }else
      {  msgs.addMessage(new Message(String.format("(QuacPac) 0 tautomer found: %s\n", preSmi),
            Message.Level.COMMENT, null));
      }
      tmpMol.delete();

//System.err.printf("%s>>%s\n",preSmi,postSmi);
//System.err.printf("After taut: %s\n", OETools.molToCanSmi(in, true));
      return true;
   }


   @Override
   public void delete() {
      tEvaluator.delete();
      tautomerOptions.delete();
   }
}
