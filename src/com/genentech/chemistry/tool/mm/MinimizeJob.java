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
package com.genentech.chemistry.tool.mm;

public class MinimizeJob
{
   private String moleculeName   = "";
   private String inputFilename  = "";
   private String outputFilename    = "";
   private String fixedAtomsIndices = "";
   private String torsionSMARTS     = null;
   
   MinimizeJob (String moleculeName, String inputFilename, String outputFilename, String fixedAtomIndices)
   {
      this.setMoleculeName(moleculeName);
      this.setInputFilename(inputFilename);
      this.setOutputFilename(outputFilename);
      this.setFixedAtomsIndices(fixedAtomIndices);
   }

   public String getMoleculeName()
   {
      return moleculeName;
   }

   public void setMoleculeName(String moleculName)
   {
      this.moleculeName = moleculName;
   }

   public String getInputFilename()
   {
      return inputFilename;
   }

   public void setInputFilename(String inputFilename)
   {
      this.inputFilename = inputFilename;
   }

   public String getOutputFilename()
   {
      return outputFilename;
   }

   public void setOutputFilename(String outputFilename)
   {
      this.outputFilename = outputFilename;
   }

   public String getFixedAtomsIndices()
   {
      return fixedAtomsIndices;
   }

   public void setFixedAtomsIndices(String fixedAtomsIndices)
   {
      this.fixedAtomsIndices = fixedAtomsIndices;
   }

   public void setTorsionSMARTS(String torsionSMARTS)
   {
      this.torsionSMARTS      = torsionSMARTS;      
   }
   
   public String getTorsionSMARTS()
   {
      return this.torsionSMARTS;      
   }
}
