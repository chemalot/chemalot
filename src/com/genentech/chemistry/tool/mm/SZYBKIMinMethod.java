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
import java.util.HashMap;
import java.util.Map;

import com.genentech.oechem.tools.OETools;

import openeye.oechem.OEAtomBase;
import openeye.oechem.OEAtomBaseIter;
import openeye.oechem.OEGraphMol;
import openeye.oechem.oechem;
import openeye.oeszybki.*;
import openeye.oechem.oemolithread;
import openeye.oechem.oemolothread;

/*
 * Implementation of a Minimization Method using the SZYBKI program from Openeye
 */
public class SZYBKIMinMethod extends MMMinMethod
{
   private final String    NAME       = "SZYBKI";
   private final double    FORCE_CONSTRAINT = 15700.0;

   @Override
   public String getMethodName()
   {
      return this.NAME;
   }

   @Override
   protected  void loadForcefielNameCmdLineMap()
   {
      addAvailableForceFieldName ("MMFF94",  null);
      addAvailableForceFieldName ("MMFF94S", "-MMFF94S");
   }

   @Override
   protected  void loadSolventNameCmdLineMap()
   {
      addAvailableSolventName ("VACUUM", null);
      addAvailableSolventName("SHEFFIELD", "-sheffield");
   }

   /*
    * Constructor
    */
   SZYBKIMinMethod ()
   {
      this.setAddIsotopeProperty(true);
   }


   @Override
   public void executeJob(MinimizeJob job, File workDir)
   {
      String tempOutFilename = job.getInputFilename().substring(0, job.getInputFilename().lastIndexOf(".")) + "_out.sdf";
      job.setOutputFilename(tempOutFilename);

      // SZYBKI options
      OESzybki         sz = new OESzybki();
      OESzybkiResults res = new OESzybkiResults();

      setSZYBKIForceFieldOptions(sz);

      // Start input stream on inputfile
      oemolithread iMolThread = new oemolithread(job.getInputFilename());
      OEGraphMol mol = new OEGraphMol();

      // Start output stream on outputfile
      oemolothread oMolThread = new oemolothread(job.getOutputFilename());

      // For molecules in inputfile
      while (oechem.OEReadMolecule(iMolThread, mol))
      {
         // Save old SD data for some reason

         OEGraphMol oldMol = new OEGraphMol();
         oechem.OECopySDData(oldMol, mol);

         setConstraints(mol, job, sz);

         // Minimize
         sz.call(mol, res);

         // Write minimized mol
         oechem.OECopySDData(mol, oldMol);

         String totalEnergyStr = String.format("%4.4f", res.GetTotalEnergy());

         oechem.OESetSDData(mol, "Total_energy", totalEnergyStr);

         // Cleanup
         removeIsotopeProperties(mol);
         removeMapIndices(mol);
         
         oechem.OEWriteMolecule(oMolThread, mol);     

         res.Print(oechem.getOeerr());
      }

      oMolThread.close();
      iMolThread.close();

      res.delete();
      sz.delete();

   }

   /*
    * Remove the isotope properties
    */
   private static void removeIsotopeProperties(OEGraphMol mol)
   {
      for (OEAtomBase at : mol.GetAtoms())
      {
         at.SetIsotope(0);
      }
   }
   
   /*
    * Remove the map index properties
    */
   private static void removeMapIndices(OEGraphMol mol)
   {
      for (OEAtomBase at : mol.GetAtoms())
      {
         at.SetMapIdx(0);
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

   /*
    * Set the torsion constraints on SZYBKI
    * @param mol
    * @param job
    * @param sz
    */
   private void setConstraints(OEGraphMol mol, MinimizeJob job, OESzybki sz)
   {
      // Get fixed atom tag value
      String atomTagValue = job.getFixedAtomsIndices();

      // Check for atoms to be fixed
      OEAtomBase[] fixAtoms = null;

      if (atomTagValue != null && atomTagValue.length()>0)
      {
         // Get the atoms from the job definition
         String[]   strIndices = atomTagValue.split(" ");
         int[]    nAtomIndices = getIntsFromStrings(strIndices);

         fixAtoms = getAtomsFromAtomIndices(mol, nAtomIndices);
      }

      if (fixTorsion() && fixAtoms != null && fixAtoms.length == 4)
      {
         // Constrain to the current torsion angle
         OEAtomBase [] sortedAtoms = sortTorsionAtomsByConnectivity(fixAtoms);
         double rads = oechem.OEGetTorsion(mol, sortedAtoms[0], sortedAtoms[1], sortedAtoms[2], sortedAtoms[3]);

         String msg = String.format("\nFixing torsion: atoms with atomic numbers: %s %s %s %s at angle: %f\n",
                                       sortedAtoms[0].GetAtomicNum(), sortedAtoms[1].GetAtomicNum(), sortedAtoms[2].GetAtomicNum(), sortedAtoms[3].GetAtomicNum(),  rads * (360.0 / (2*3.1415)) );
         System.err.print(msg);
         System.err.print("SMARTS: " + job.getTorsionSMARTS() + "\n");

         sz.SetTorsionConstraint(job.getTorsionSMARTS(), FORCE_CONSTRAINT, rads);
      }
   }


   private void setSZYBKIForceFieldOptions(OESzybki sz)
   {
      // Forcefield
      int szybkiFF = OEForceFieldType.MMFF94;
      String ffChoice = getForceFieldChoice();
      if (ffChoice != null)
      {
         if (ffChoice.equalsIgnoreCase("MMFF94S"))
         {
            szybkiFF = OEForceFieldType.MMFF94;
         }
      }
      sz.SetForceFieldType(szybkiFF);

      // Solvent
      int szybkiSolvent = OESolventModel.NoSolv;
      String solvChoice = getSolventChoice();
      if (solvChoice != null)
      {
         if (solvChoice.equalsIgnoreCase("SHEFFIELD"))
         {
            szybkiSolvent = OESolventModel.Sheffield;
         }
      }
      sz.SetSolventModel(szybkiSolvent);
   }

}
