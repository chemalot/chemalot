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
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import openeye.oechem.OEAtomBase;
import openeye.oechem.OEAtomBaseIter;
import openeye.oechem.OEGraphMol;
import openeye.oechem.OEMolBase;
import openeye.oechem.OESMILESFlag;
import openeye.oechem.OEUnaryAtomBoolFunc;
import openeye.oechem.OEUnaryAtomPred;
import openeye.oechem.oechem;
import openeye.oechem.oemolithread;
import openeye.oechem.oemolothread;

import org.apache.cxf.common.util.StringUtils;

import com.aestel.utility.SimpleProcess;

public abstract class MMMinMethod
{
   protected static final String PROGRAM_TAG    = "Program";
   protected static final String FORCEFIELD_TAG = "ForceField";
   protected static final String METHOD_TAG     = "MinMethod";
   protected static final String SOLVENT_TAG    = "Solvent";
  
   private String forceFieldChoice  = "";
   private String solventChoice     = "";
   
   private final Map<String,String> forceFieldNameCmdLineMap = new HashMap<String,String>();
   private final Map<String,String> solventNameCmdLineMap    = new HashMap<String,String>();  
   private boolean addIsotopeProperty = false;
   private boolean fixTorsion = false;
   
   private static class PredAtomHasIsotope extends OEUnaryAtomPred 
   {
      public boolean constCall(OEAtomBase atom) 
      {
          if (atom.GetIsotope() != 0)
          {
             return true;
          }
          return false;
      }
      public OEUnaryAtomBoolFunc CreateCopy() {
          OEUnaryAtomBoolFunc copy = new PredAtomHasIsotope();
          copy.swigReleaseOwnership();
          return copy;
      }
  }

   
   /*
    * Default constructor
    * Load the available forcefield names and solvents
    */
   public MMMinMethod()
   {
      loadForcefielNameCmdLineMap();
      loadSolventNameCmdLineMap();
   }
   
   /*
    * Required method to load the available forcefield names 
    *  explicit command line options for those forcefields
    * @return void
    */
   protected abstract void loadForcefielNameCmdLineMap();
   
   /*
    * Required method to load the available solvent names and the 
    *  explicit command line options for those forcefields
    *  @return void
    */
   protected abstract void loadSolventNameCmdLineMap();
  
   /*
    * Return the base name of the minimize program (e.g. SZYBKI)
    * @return the String name of the program used
    */
   public abstract String getMethodName();
   
   /*
    * Return an array of strings containing the list of allowed force fields
    * @return array of string forcefield names equivalent to the command line choices
    */
   public String[] getAvailableForceFieldNames()
   {
      Set<String> keys = forceFieldNameCmdLineMap.keySet();
      return keys.toArray(new String[keys.size()]);
   }
   
   /*
    * Return an array of strings containing the list of allowed solvents
    * @return array of string solvent names equivalent to the command line choices
    */
   public String[] getAvailableSolventNames()
   {
      Set<String> keys = solventNameCmdLineMap.keySet();
      return keys.toArray(new String[keys.size()]);
   }
   
   /*
    * Get the force field choice.
    * @return the name of the single chose forcefield  
    */
   public String getForceFieldChoice()
   {
      return this.forceFieldChoice;
   }
   
   /*
    * Get the force field command line option.  
    * @param forceFieldName the name of the forcefield
    * @return the explicit command line option for this forcefield
    * to be used internally when calling the program.
    */
   public String getForceFieldCmdLine(String forceFieldName)
   {
      if (forceFieldNameCmdLineMap.containsKey(forceFieldName) && 
          forceFieldNameCmdLineMap.get(forceFieldName) != null)
      {
         return forceFieldNameCmdLineMap.get(forceFieldName);
      }
      return null;
   }
   
   /*
    * Get the chosen solvent name.
    * @return the name of the chosen solvent  
    */
   public String getSolventChoice()
   {
      return this.solventChoice;
   }
   
   /*
    * Get the solvent command line option.
    * @param solventName the name of the solvent   
    * @return the explicit command line option for this solvent
    * to be used internally when calling the program.
    */
   public String getSolventCmdLine(String solventName)
   {
      if (solventNameCmdLineMap.containsKey(solventName) && 
          solventNameCmdLineMap.get(solventName) != null)
      {
         return solventNameCmdLineMap.get(solventName);
      }
      return null;
   }
   
   /*
    * Get flag 
    * @return true if isotope property data should be added to torsion atoms.
    */
   public boolean addIsotopeProperty()
   {
      return addIsotopeProperty;
   }

   /*
    * Set whether isotope property data should be added to torsion atoms
    */
   public void setAddIsotopeProperty(boolean addIsotopeProperty)
   {
      this.addIsotopeProperty = addIsotopeProperty;
   }

   /*
    * Execute the minimization method
    * @param inFilename the input filename to be minimized
    * @param outFilename the output filename to which minimized results will be written
    * @param fixedAtomTag the SDF tag in inFilename which holds the indices of for atoms to be held fixed
    * @param workDirPath the working directory path
    * @return void
    */
   public void execute(String inFilename, String outFilename, String fixedAtomTag, boolean fixTorsions, String workDirPath)
   {
      // Check working directory
      File workDir = new File (workDirPath);
      if (!workDir.exists())
      {
         throw new Error("Cannot find working directory " + workDirPath);
      } 
      
      setFixTorsion(fixTorsions);
      
      String prefix = "";     
      
      // Split the input sd file into separate temp files based on the fixed atom tag.
      List<MinimizeJob> jobs = 
               SDFSplitOnTagValue.splitOnTagValue(inFilename, fixedAtomTag, prefix, workDirPath);
      
      if (addIsotopeProperty)
      {
         addIsotopeProperty(jobs);
      }
      
      // For each temp file
      List<String> tempOutputFilenames = new ArrayList<String>();
      for (MinimizeJob job : jobs)
      {         
         // Get command line         
         executeJob (job, workDir);
                  
         // Save the output file
         tempOutputFilenames.add(job.getOutputFilename());
      }
      
      Map<String,String> sdTagsToAdd = new HashMap<String, String>();
      
      sdTagsToAdd.put(PROGRAM_TAG,    this.getMethodName());
      sdTagsToAdd.put(FORCEFIELD_TAG, this.getForceFieldChoice());
      sdTagsToAdd.put(SOLVENT_TAG,    this.getSolventChoice());
      sdTagsToAdd.put(METHOD_TAG,     this.getMethodName() + "_" + this.getForceFieldChoice() + "_" + this.getSolventChoice());
            
      // Cleanup
      SDFSplitOnTagValue.concatenateFilesAndWriteTags(tempOutputFilenames, outFilename, sdTagsToAdd);     
   }
   
   protected static void addIsotopePropertyToJob(MinimizeJob job)
   {
      if (job.getFixedAtomsIndices() == null || job.getFixedAtomsIndices().length() < 1)
      {
         return;
      }
      
      // Get input filename
      String oldInputFilename = job.getInputFilename();
      
      // Create new input filename for isotope
      String newInputFilename = oldInputFilename.replaceFirst(".sdf", "_isotope.sdf");
      
      String torsionSMARTS = "";
      // Start output thread 
      oemolothread oMolIntermedInFileThread = new oemolothread(newInputFilename);
      
      // Start input thread with infile
      oemolithread iMolThread = new oemolithread(oldInputFilename);
      OEGraphMol mol = new OEGraphMol();
      
      // Read all molecules
      while (oechem.OEReadMolecule(iMolThread, mol))
      {     
         // Get the atoms from the job definition
         String[]   strIndices = job.getFixedAtomsIndices().split(" ");        
         OEAtomBase[] torAtoms = getAtomsFromAtomIndices(mol, getIntsFromStrings(strIndices));
         OEAtomBase[] sortedAtoms = sortTorsionAtomsByConnectivity (torAtoms);   
         
         // Label torsion atoms with isotopes to identify them
         int atomIdx = 1;
         for (OEAtomBase atom : sortedAtoms)
         {
            Integer isotope = atom.GetAtomicNum()*2;
            atom.SetIsotope(isotope);   
            atom.SetMapIdx(atomIdx++);
         }
         
         // If first time, then also create the torsion smarts
         if (job.getTorsionSMARTS() == null )
         {
            OEGraphMol tempTorMol = new OEGraphMol();
            oechem.OESubsetMol(tempTorMol, mol, new PredAtomHasIsotope());
            //torsionSMARTS = oechem.OECreateSmiString(tempTorMol, OESMILESFlag.Isotopes|OESMILESFlag.AtomMaps);
            torsionSMARTS = getTorsionSMARTSWithIsotopes(tempTorMol);
            job.setTorsionSMARTS(torsionSMARTS);
         }
         
         oechem.OEWriteMolecule(oMolIntermedInFileThread, mol);
      }
                        
      // Close current threads
      iMolThread.close();
      iMolThread.delete();
      oMolIntermedInFileThread.close();
      oMolIntermedInFileThread.delete();
      
      // Set new input filename
      job.setInputFilename(newInputFilename);
   }
   
   
   
   /**
    * Get the smarts for the input mol
    * @param mol
    * @return a smarts string
    */
   protected static String getTorsionSMARTSWithIsotopes(OEMolBase mol)
   {
      
      OEAtomBase [] atomList = new OEAtomBase[mol.NumAtoms()];
      String smarts = "";
      int nIndex = 0;
      for (OEAtomBase atom : mol.GetAtoms()) 
      {
         atomList[nIndex] = atom;
         nIndex++;
      }
      OEAtomBase[] sortedAtoms = sortTorsionAtomsByConnectivity (atomList);   
      
      //  Manually construct the SMARTS since simply getting the SMILES doesn't seem
      //   to preserve the atom order.
      nIndex = 1;
      for (OEAtomBase sortedAtom : sortedAtoms) 
      {         
         if (smarts.length() > 0)
         {
            smarts = smarts + "~";
         }
         smarts += String.format("[%d*:%d]",sortedAtom.GetIsotope(), nIndex); 
         nIndex++;
      }
      return smarts;
   }

   
   /**
    * Sort the incoming array of atoms by their connectivity.  First terminal atom in incoming array will be "first" in output
    * Assumes these are a torsion.
    * @param torAtoms atoms to sort
    * @return an array of atoms that are sorted by connectivity
    */
   protected static OEAtomBase[] sortTorsionAtomsByConnectivity (OEAtomBase[] torAtoms)
   {
      OEAtomBase[] sortedAtoms = new OEAtomBase[torAtoms.length];
      
      // Create a map for easy lookup
      Map<OEAtomBase,Integer> atomNeighborCountMap = new HashMap<OEAtomBase,Integer>();      
      for (OEAtomBase atom : torAtoms)
      {  atomNeighborCountMap.put(atom, 0);  }
      
      
      // Find any terminal end
      OEAtomBase currAtom = null;      
      for (OEAtomBase atom : torAtoms)
      {
         currAtom = atom;
         
         int nNbrsInTorsion = 0;         
         for (OEAtomBase nbrAtom : atom.GetAtoms())
         {
            //System.err.println(nbrAtom.GetAtomicNum());
            if (atomNeighborCountMap.containsKey(nbrAtom))
            {
               nNbrsInTorsion++;
            }
         }
         
         if (nNbrsInTorsion == 1)
         {
            break;
         }
      }
      
      // Now start at terminal atom and save atoms in order of torsion connectivity
      int nSortIndex = 0;
      int failsafe   = 0;
      sortedAtoms[nSortIndex] = currAtom;
      nSortIndex++;
      boolean foundFour = false;
      OEAtomBase prevAtom = null;
      
      while (nSortIndex<4 && failsafe<100)
      {            
         // Find neighbor in torsion
         for (OEAtomBase nbrAtom : currAtom.GetAtoms())
         {
            //System.err.println(nbrAtom.GetAtomicNum());
            if (   atomNeighborCountMap.containsKey(nbrAtom)  
              && ((prevAtom == null) || ( prevAtom != null && (nbrAtom.GetIdx() != prevAtom.GetIdx())))
               )
            {
               // Got it
               sortedAtoms[nSortIndex] = nbrAtom;
               nSortIndex++;
               
               prevAtom = currAtom;
               currAtom = nbrAtom;
               break;
            }           
         }
                 
         failsafe++;         
      }
      
      if (failsafe==100)
      {
         throw new Error("Couldn't find torsion atoms when sorting by connectivity.\n");
      }
      
      return sortedAtoms;
   }

   
   /**
    * Add the isotope property to the atoms which are specified in the 
    * input file in the jobs objects.  Write new input files with these
    * isotope properties and set them as the new input files.
    * @param jobs the job instances which contain the file information
    */
   protected static void addIsotopeProperty(List<MinimizeJob> jobs)
   {
      for (MinimizeJob job : jobs)
      {
         addIsotopePropertyToJob(job);
      }
   }

   
   /**
    * Internal method to get atom object references from an array of atom indices into a molecule
    * @param mol the molecule to search for atoms
    * @param nAtomIndices the array of zero-based int atom indices
    * @return an array of OEAtomBase object references
    */
   protected static OEAtomBase[] getAtomsFromAtomIndices(OEGraphMol mol, int[] nAtomIndices)
   {
      OEAtomBase[] atoms = new OEAtomBase[nAtomIndices.length];
      int returnArrIndex = 0;
      for (int nAtomIndex : nAtomIndices)
      {
         // zero-based index added to iterator starting at first molecule
         OEAtomBaseIter aItr = mol.GetAtoms();
         OEAtomBase     atom = aItr.ToFirst().Increment(nAtomIndex).Target();   
         atoms[returnArrIndex++] = atom;
      }                                        
      return atoms;   
   }

   /**
    * Get an array if ints from an array of strings
    * @param strIndices
    * @return
    */
   protected static int[] getIntsFromStrings(String[] strIndices)
   {
      int[] retInts = new int[strIndices.length];
      int returnArrIndex = 0;
      for (String strIndex : strIndices)
      {
         retInts[returnArrIndex++] = Integer.parseInt(strIndex);
      }
      return retInts;
   }

   /*
    * Execute the job for the specific implementation
    * @return the command line string to be executed for this method
    */
   protected abstract void executeJob(MinimizeJob job, File workDir);
   
   /*
    * Store an available forcefield name and the command line option equivalent
    * @param name the forcefield name used 
    * @param cmdLine the explicit option to be use when calling the internal program
    * @return void
    */
   protected final void addAvailableForceFieldName(String name, String cmdLine)
   {
      this.forceFieldNameCmdLineMap.put(name, cmdLine);
   }
   
   /*
    * Set an available solvent name ad the cmdline option equivalent
    * @param name the solvent name used 
    * @param cmdLine the explicit option to be use when calling the internal program
    * @return void
    */
   protected final void addAvailableSolventName(String name, String cmdLine)
   {
      this.solventNameCmdLineMap.put(name, cmdLine);
   }
   
   /*
    * Internal method to fire off command on command line cshell commandline
    * @param commnd the method to call
    * @return void
    */
   protected static void executeCommand (String command)
   {
      System.err.printf("Executing: " + command.toString());      
     
      try
      {
         SimpleProcess p = new SimpleProcess(true, null, "csh", "-c", command.toString());
         if (!StringUtils.isEmpty(p.getProgramOutput()))
            System.err.println(p.getProgramOutput());
         if (!StringUtils.isEmpty(p.getProgramErrorOutput()))
            System.err.println(p.getProgramErrorOutput());
      }
      catch(Exception e)
      {
         throw new Error(e);
      }
   }
     

   /*
    * Set the forcefield to use for this execution
    * @param requestedForceField the name of the requested forcefield
    * @return void
    */
   public void setForceField(String requestedForceField)
   {      
      boolean foundNewOption = false;
      for (String availableFF : getAvailableForceFieldNames())
      {
         if (availableFF.equalsIgnoreCase(requestedForceField.toLowerCase()))
         {
            this.forceFieldChoice = requestedForceField;
            foundNewOption = true;
            break;
         }
      }
      if (!foundNewOption)
      {
         throw new Error ("Forcefield choice " + requestedForceField + " not available for program " + getMethodName());
      }      
   }
   
   /*
    * Set the solvent to use for this execution
    * @param requestedSolvent the name of the requested forcefield
    * @return void
    */
   public void setSolvent(String requestedSolvent)
   {
      boolean foundNewOption = false;
      for (String availableSolvent : getAvailableSolventNames())
      {
         if (availableSolvent.equalsIgnoreCase(requestedSolvent.toLowerCase()))
         {
            this.solventChoice = requestedSolvent;
            foundNewOption = true;
            break;
         }
      }
      if (!foundNewOption)
      {
         throw new Error ("Solvent choice " + requestedSolvent + " not available for program " + getMethodName());
      }      
   }

   public boolean fixTorsion()
   {
      return fixTorsion;
   }

   public void setFixTorsion(boolean fixTorsion)
   {
      this.fixTorsion = fixTorsion;
   }

}
