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
package com.genentech.retrival.SDFExport;

import java.util.regex.Pattern;

public class SDRecord {

   private final StringBuilder data;
   private final String mol;
   private static final Pattern NLEndPat = Pattern.compile("[\r\n]+$");
   private static final Pattern REMOVEDollars = Pattern.compile("\\s+\\$\\$\\$\\$\\s*$");

   static final String NEWLINE = "\r\n";  // windows stile for now unix: "\n"

   public SDRecord(String mol, String data) {
      // normalize newline characters to unix mode
      mol  = mol.replace("\r\n", "\n").replace("\n\r", "\n");
      data = data.replace("\r\n", "\n").replace("\n\r", "\n");

      mol = REMOVEDollars.matcher(mol).replaceAll("");

      //normalize newline characters to windows.
      mol = mol.replace("\n", NEWLINE);
      data = data.replace("\n", NEWLINE);
      this.mol = NLEndPat.matcher(mol).find() ? mol : mol + NEWLINE; // ensure newline at end
      this.data = new StringBuilder(data);
   }

   /**
    * all data fields in sdf format.
    * not including the final $$$$.
    */
   public String getData() {
      return data.toString();
   }

   public void addData(String tag, String value) {
      data.append("> <").append(tag).append(">").append(NEWLINE)
          .append(value).append(NEWLINE);
      if( value != null && value.length() > 0 )
         data.append(NEWLINE);
   }

   /**
    * Get molfile.
    */
   public String getMolStr() {
      return mol;
   }

   /**
    * sdf record without terminating "\n".
    */
   @Override
   public String toString() {
      return mol+data+"$$$$";
   }
}
