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

public class MOEMinMethod extends MMMinMethod
{   
   private final String    NAME       = "MOE";
   
   @Override
   public String getMethodName()
   {
      return this.NAME;
   }
   
   @Override
   protected  void loadForcefielNameCmdLineMap()
   {   
      addAvailableForceFieldName ("AMBER12EHT", "Amber12EHT");
      addAvailableForceFieldName ("CHARM22",    "charm22");
      addAvailableForceFieldName ("CHARM27",    "charm27");
      addAvailableForceFieldName ("EMPIRICAL",  "empirical");
      addAvailableForceFieldName ("ENGH_HUBER", "engh_huber");
      addAvailableForceFieldName ("KOLL89",     "koll89");
      addAvailableForceFieldName ("KOLL94",     "koll94");
      addAvailableForceFieldName ("KOLL99",     "koll99");
      addAvailableForceFieldName ("MMFF94",     "mmff94");
      addAvailableForceFieldName ("MMFF94S",    "mmff94s");
      addAvailableForceFieldName ("MMFF94X",    "mmff94x");
      addAvailableForceFieldName ("OPLSAA",     "oplsaa");
      addAvailableForceFieldName ("PEF95sac",   "pef95sac");
      addAvailableForceFieldName ("PFROSST",    "pfrosst");
      addAvailableForceFieldName ("TAFF",       "taff");
   }
 
   @Override
   protected  void loadSolventNameCmdLineMap()
   {
      addAvailableSolventName ("VACUUM",   "Vacuum");
      addAvailableSolventName ("DISTANCE", "Distance");
      addAvailableSolventName ("R-Field",  "R-Field");   
      addAvailableSolventName ("BORN",     "Born");
   }
      
   @Override
   protected void executeJob (MinimizeJob job, File workDir)
   {        
      int nBeginFilenameIndex = 0;
      if (job.getInputFilename().lastIndexOf("/") > -1)
      {
         nBeginFilenameIndex = job.getInputFilename().lastIndexOf("/")+1;
      }
      else if (job.getInputFilename().lastIndexOf("\\") > -1)
      {
         nBeginFilenameIndex = job.getInputFilename().lastIndexOf("\\")+1;
      }
      String tempOutFilename = job.getInputFilename().substring(nBeginFilenameIndex, job.getInputFilename().lastIndexOf(".")) + "_out.sdf";
      String workPath =  workDir.getAbsolutePath();
      tempOutFilename = workPath + "/" + tempOutFilename;
      job.setOutputFilename(tempOutFilename);
           
      // Build up command line
      StringBuilder command = new StringBuilder();
      command.append("cd " + workDir.getAbsolutePath() + ";");    
      
      // Temporary beta version if MOE which has some fixes for sdminimize
      command.append("setenv PATH /gne/home/benjamds/workspace/programs/moe2014.0901/bin:$PATH" + ";");
      command.append("setenv MOE /gne/home/benjamds/workspace/programs/moe2014.0901" + ";");
      command.append("cat " + job.getInputFilename() + " | sdminimize -o " + tempOutFilename);                
      
      String ffCmdLine = getForceFieldCmdLine(getForceFieldChoice());
      if (ffCmdLine != null)
      {
         command.append(" -ff " + ffCmdLine);
      }
      
      String solvCmdLine = getSolventCmdLine(getSolventChoice());
      if (solvCmdLine != null)
      {
         command.append(" -solvation " + solvCmdLine);
      }
      
      if (job.getFixedAtomsIndices() != null)
      {
         String moeFixedAtomIndices = convertIndicesStringToMOEString(job.getFixedAtomsIndices());
         
         if (fixTorsion())
         {
            command.append(" -restrainDihAtoms " + moeFixedAtomIndices);
            command.append(" -restrainforce 5000");
         }
         else
         {
            command.append(" -fixAtoms " + moeFixedAtomIndices);
         }
      }

      // Execute
      executeCommand(command.toString());  
   }
    
   /*
    * Internal helper to convert the zero-based string of atom indices coming from the input 
    * sdf file to a one's-based string of indices as needed by MOE.
    * @param zeroBasedIndexString the zero-based index of atoms "0 1 2 4"
    * @return a string of one's based indices "1 2 3 5"
    */
   private static String convertIndicesStringToMOEString(String zeroBasedIndexString)
   {
      String [] indexStrings = zeroBasedIndexString.split("\\s+");
      String newIndexString = "";
      for (String indexStr : indexStrings)
      {
         Integer index = Integer.valueOf(indexStr);
         
         // MOE wants indices starting at 1 but we have 0-based index
         index++;
         
         newIndexString += index.toString() + ",";        
      }
      newIndexString = newIndexString.substring(0,newIndexString.lastIndexOf(","));
      return newIndexString;
   }
   
}
