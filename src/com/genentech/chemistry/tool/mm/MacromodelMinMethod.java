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

import java.io.File;

/*
 * Class implementation representing a minimization method
 * using the Schrodinger program, Macromodel
 */
public class MacromodelMinMethod extends MMMinMethod
{
   private final String    NAME       = "MACROMODEL";       
      
   @Override
   public String getMethodName()
   {
      return this.NAME;
   }

   @Override
   protected  void loadSolventNameCmdLineMap()
   {
      // Solvent methods
      addAvailableSolventName("VACUUM",    "None");
      addAvailableSolventName("WATER",     "Water");
      addAvailableSolventName("OCTANOL",   "Octanol");   
      addAvailableSolventName("CHCL3",     "CHCL3");
   }
   
   @Override
   protected  void loadForcefielNameCmdLineMap()
   {
   // Forcefields
      addAvailableForceFieldName("MM2",       "MM2*");
      addAvailableForceFieldName("MM3",       "MM3*");
      addAvailableForceFieldName("AMBER",     "AMBER*");
      addAvailableForceFieldName("AMBER94",   "AMBER94");
      addAvailableForceFieldName("OPLS",      "OPLS");
      addAvailableForceFieldName("MMFF",      "MMFF");
      addAvailableForceFieldName("MMFFS",     "MMFFs");
      addAvailableForceFieldName("OPLS_2001", "OPLS_2001");
      addAvailableForceFieldName("OPLS_2005", "OPLS_2005");
      addAvailableForceFieldName("OPLS_21",   "OPLS_21");
   }
 
 
   @Override
   protected void executeJob (MinimizeJob job, File workDir)
   {      
      String tempOutFilename = job.getInputFilename().substring(0, job.getInputFilename().lastIndexOf(".")) + "_out.sdf";
      job.setOutputFilename(tempOutFilename);
           
      // Build up command line
      StringBuilder command = new StringBuilder();
      command.append("cd " + workDir.getAbsolutePath() + ";");       
      command.append("$SCHRODINGER/run macromodelmin.py " + job.getInputFilename() + " " + tempOutFilename);                     
      
      String ffCmdLine = getForceFieldCmdLine(getForceFieldChoice());
      if (ffCmdLine != null)
      {
         command.append(" -forcefield " + ffCmdLine);
      }
      
      String solvCmdLine = getSolventCmdLine(getSolventChoice());
      if (solvCmdLine != null)
      {
         command.append(" -solvent " + solvCmdLine);
      }
      
      if (job.getFixedAtomsIndices() != null)
      {
         String macromodelFixedAtomIndices = convertIndicesStringToMacromodelString(job.getFixedAtomsIndices());
         
         if (fixTorsion())
         {
            command.append(" -constrain_by_index " + macromodelFixedAtomIndices);
         }
         
      }
                          
      // Execute
      executeCommand(command.toString());      
   }
 
   
   private static String convertIndicesStringToMacromodelString(String zeroBasedIndexString)
   {
      String [] indexStrings = zeroBasedIndexString.split("\\s+");
      String newIndexString = "";
      for (String indexStr : indexStrings)
      {                
         Integer index = Integer.valueOf(indexStr);
         
         // MACROMODEL wants indices starting at 1 but we have 0-based index
         index++;
         
         newIndexString += index.toString() + ",";    
      }
      // remove last comma
      newIndexString = newIndexString.substring(0,newIndexString.lastIndexOf(","));
      return newIndexString;
   }

 
}
