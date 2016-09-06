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

import openeye.oechem.OEGraphMol;
import openeye.oechem.oechem;

import org.jdom.Element;

import com.aestel.utility.MessageList;


/**
 * Abstract class implementing basic functionality for structure checking rules.
 * 
 * @author albertgo
 *
 */
public abstract class AbstractStructureCheck implements StructureCheckInterface{

   protected final String msgTxt;
   private final String description;
   private final String example;
   private final HydrogenMode reqHydrogenMode;
   private final String checkName;

   public AbstractStructureCheck(Element elem) {
      msgTxt = elem.getChildTextTrim("message");
      description = elem.getChildTextTrim("description");
      example = elem.getChildTextTrim("example");
      checkName = elem.getName();
      
      String hModeStr = elem.getAttributeValue("reqHydrogenMode");
      if(hModeStr != null) {
         reqHydrogenMode = HydrogenMode.valueOf(hModeStr);
      }else {
         reqHydrogenMode = HydrogenMode.ANY;
      }
   }
   
   public AbstractStructureCheck(Element elem, HydrogenMode reqHydrogenMode) {
      super();

      msgTxt = elem.getChildTextTrim("message");
      description = elem.getChildTextTrim("description");
      example = elem.getChildTextTrim("example");
      checkName = elem.getName();
      this.reqHydrogenMode = reqHydrogenMode;
   }

   @Override
   public String getDescription() {
      return description;
   }

   @Override
   public String getExampleInput() {
      return example;
   }

      
   /**
    * Default implementation of {@link StructureCheckInterface#checkExample()}
    * simply call {@link #checkExample()} with the example structure and return
    * true if at least one message was created.
    * 
    * This may be overwritten to implement more specific checking.
    */
   @Override
   public boolean checkExample() {
      // no example defined
      if(getExampleInput() == null) return true;
      
      OEGraphMol mol = new OEGraphMol();
      oechem.OEParseSmiles(mol, getExampleInput());
      
      MessageList msgs = new MessageList();
      checkStructure(mol, null, msgs );
      mol.delete();
      

      // check that there is something to report
      return msgs.countMessages() > 0;
   }
   
   @Override
   public HydrogenMode getRequiredHydrogenMode() {
      return reqHydrogenMode;
   }
   
   @Override
   public String getCheckName() {
      return checkName;
   }
   @Override
   public String getDescriptionHTML() {
      return getDescription();
   }

   // TODO this should be elsewhere
   public static String encodeHTML(String s) {
       StringBuilder out = new StringBuilder();
       for(int i=0; i<s.length(); i++) {
           char c = s.charAt(i);
           if(c > 127 || c=='"' || c=='<' || c=='>')
              out.append("&#"+(int)c+";");
           else
               out.append(c);
       }
       return out.toString();
   }
}
