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
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.UUID;

import openeye.oechem.OEGraphMol;
import openeye.oechem.oechem;
import openeye.oechem.oemolithread;
import openeye.oechem.oemolothread;

/*
 * Split the input file into multiple temp files grouped by the tag value
 * TODO: First implementation assumes tag values are pre-grouped together.  Make generic 
 */
public class SDFSplitOnTagValue
{
   /*
    * @return a List of job instances which define the input/output filenames
    */
   public static List<MinimizeJob> splitOnTagValue (String inputSDFile, String tagName, String prefix, String workDirPath)
   {
      List<MinimizeJob> jobs = new ArrayList<MinimizeJob>();
      
      // Check working directory
      File workDir = new File (workDirPath);
      if (!workDir.exists())
      {
         throw new Error("Cannot find working directory " + workDirPath);
      }           
      
      String   runPrefix      = prefix + UUID.randomUUID().toString();
      Integer  nIntermedIndex = 0;
      
      String   outputTempFilePath = workDirPath + "/" + runPrefix + "_input_" + nIntermedIndex.toString() + ".sdf";   
      
      String   cachedTagValue = null;      
      String   cachedMolName  = null;
      
      // Start output thread with temp filename
      oemolothread oMolIntermedInFileThread = new oemolothread(outputTempFilePath);
      
      // Start input thread with infile
      oemolithread iMolThread = new oemolithread(inputSDFile);
      OEGraphMol mol = new OEGraphMol();
      
      // Read all molecules
      while (oechem.OEReadMolecule(iMolThread, mol))
      {     
         // Get sdf tag value
         String currentTagValue = null;
         
         if (tagName != null && oechem.OEHasSDData(mol, tagName))
         {
            currentTagValue = oechem.OEGetSDData(mol, tagName);
         }
         
         // if sdf tag is different than last, 
         if ( currentTagValue  != null &&
              cachedTagValue   != null && 
             !currentTagValue.trim().equalsIgnoreCase(cachedTagValue.trim())
            )
         {
            //
            // We have a new tag value string which means a different temp file.
            //   so close down the current thread, save the previous file and start a new 
            //   output temp sdf file.
            //
            
            // Close current output thread
            oMolIntermedInFileThread.close();
            oMolIntermedInFileThread.delete();                      
        
            jobs.add(new MinimizeJob(cachedMolName, outputTempFilePath, null, cachedTagValue));            
            
            // Open new output thread
            nIntermedIndex++;
            outputTempFilePath = runPrefix + "_input_" + nIntermedIndex.toString() + ".sdf";                       
            oMolIntermedInFileThread = new oemolothread(outputTempFilePath);
      
            // Write out molecule to othread            
            oechem.OEWriteMolecule(oMolIntermedInFileThread, mol);
            
            cachedMolName = mol.GetTitle();
            cachedTagValue = currentTagValue; 
         }
         else
         {            
            // Write out molecule to othread            
            oechem.OEWriteMolecule(oMolIntermedInFileThread, mol);
            
            if (cachedTagValue == null)
            {  
               // First time
               cachedTagValue = currentTagValue; 
               cachedMolName  = mol.GetTitle();
            }
         }
      }
      
      
      // Close current threads
      iMolThread.close();
      iMolThread.delete();
      oMolIntermedInFileThread.close();
      oMolIntermedInFileThread.delete();
      
      // Save last file      
      jobs.add(new MinimizeJob(cachedMolName, outputTempFilePath, null, cachedTagValue));
      return jobs;
   }
  
   /*
    * Helper function to combine a set of files into one SDF 
    * Also add any additional tags.
    */
   static void concatenateFilesAndWriteTags(List<String> tempOutputFilenames, String outputFilename, Map<String,String> tagValueMap)
   {            
      // Concatenate the output files and pass to the outfile
      oemolothread oMolOutileThread = new oemolothread(outputFilename);
      for (String inputFilename : tempOutputFilenames)
      {
      // Start input thread with infile
         oemolithread iMolTempFileThread = new oemolithread(inputFilename);
         OEGraphMol tempMol = new OEGraphMol();
         // Read all molecules
         while (oechem.OEReadMolecule(iMolTempFileThread, tempMol))
         {  
            for (String tagName : tagValueMap.keySet())
            {
               oechem.OESetSDData( tempMol, tagName, tagValueMap.get(tagName));
            }
            
            oechem.OEWriteMolecule(oMolOutileThread, tempMol);
         }
         iMolTempFileThread.close();
         iMolTempFileThread.delete();
      }   
      oMolOutileThread.close();
      oMolOutileThread.delete();
   }
   
}
