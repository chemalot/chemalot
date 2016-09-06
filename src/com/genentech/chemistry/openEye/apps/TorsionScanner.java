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
package com.genentech.chemistry.openEye.apps;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import openeye.oechem.*;
import openeye.oeszybki.*;

import com.genentech.chemistry.openEye.conformerSampler.ConformerSampler;
import com.genentech.chemistry.openEye.util.AtomListFunctor;
import com.genentech.oechem.tools.AtomIntPropertyFunctor;
import com.genentech.oechem.tools.OETools;


/**
 * @author Ben Sellers / Alberto Gobbi / 2014
 * Copyright 2014 Genentech
 *
 * Class to support generation of torsion conformers.
 * The scanner will generate conformers and angular increments
 * and can also generate alternative conformers for atoms outside the
 * torsion bond, can perform minimizations of these atoms or do both or none.
 */
public class TorsionScanner
{
   @SuppressWarnings("unused")
   private static final String MY_NAME = "TorsionScanner";
   public  static final int FRAGNumTag = oechem.OEGetTag("fragNum");

   private static final String TORSION_SCANNER_ANGLE_FILENAME = "FixedTorsionExpansionPatterns.txt";
   public static final String ANGLE_TAG                      = "ScanVar_1";

   private final String  torsionAtomsTag;
   private final String  coreFilename;
   private final int     nSteps;
   private final double  nStartTor;
   private final double  nTorIncr;
   private final int     nMaxConfsPerStep;
   private final boolean doMinimize;
   private final String  constraintStrength;

   private final OEMolBase      bondMol;
   private final OEAtomBase[] bondAtoms;

   private ConformerSampler conformerSampler;

   private final OESzybki szybki;

   private List<OEGraphMol> cachedInputMols =  new ArrayList<OEGraphMol>();


   /**
    * Constructor
    * @param bondFile sdf file with four atoms whose positions define the bond in the input file
    * @param coreFilename core filename to write out
    * @param torsionAtomsTag the tag name to which to write out atom indices for the torsion
    * @param nSteps number of angular steps to generate around torsion
    * @param nStartTor the starting torsion angle in degrees
    * @param nTorIncr the degree increment between generated torsion angles
    * @param nMaxConfsPerStep the max number of conformers to generate using rotatable bonds other than the fixed torsion
    * @param doMinimize minimize at each step.  if the nMaxConfsPerStep is > 1, write out only the min E conformer per step
    * @param constraintStrength "strong" "medium" "weak" "none"
    */
   public TorsionScanner(
            String bondFilename,
            String coreFilename,
            String torsionAtomsTag,
            int nSteps,
            double nStartTor,
            double nTorIncr,
            int nMaxConfsPerStep,
            boolean doMinimize,
            String constraintStrength)
   {
      this.torsionAtomsTag    = torsionAtomsTag;
      this.coreFilename       = coreFilename;
      this.nSteps             = nSteps;
      this.nStartTor          = nStartTor;
      this.nTorIncr           = nTorIncr;
      this.nMaxConfsPerStep   = nMaxConfsPerStep;
      this.doMinimize         = doMinimize;
      this.constraintStrength = constraintStrength;

      // Set torsion bond mol object
      this.bondMol  = readBondAtomsFromFile(bondFilename);

      // Set torsion bond atoms
      this.bondAtoms = getOrderedBondAtoms(bondMol);

      // Create a conformerSampler for generating new conformers
      if (nMaxConfsPerStep > 1)
      {
         initializeConformerSampler();
      }

      // Initialize szybki
      if (doMinimize)
      {
         szybki = initializeSZYBKI(constraintStrength);
      }
      else
      {
         szybki = null;
      }
   }

   public String getTorsionAtomsTag()
   {  return torsionAtomsTag;  }

   public int getSteps()
   {  return nSteps;   }

   public double getTorIncr()
   {  return nTorIncr; }

   public int getMaxConfsPerStep()
   {  return nMaxConfsPerStep;   }

   public boolean doMinimize()
   {  return doMinimize;  }

   public String getConstraintStrength()
   {  return constraintStrength; }


   /*
    * Cleanup internal objects
    */
   public void close()
   {
      if (szybki != null)
      {
         szybki.delete();
      }
      if (bondMol != null)
      {
         bondMol.delete();
      }
      if (conformerSampler != null)
      {
         conformerSampler.close();
      }
      clearCore();
   }


   /**
    * Clear out the input molecule cache when computing the common core.
    */
   public void clearCore()
   {
      for( OEGraphMol m : cachedInputMols)
         m.delete();
      cachedInputMols.clear();
   }

   /**
    * for each molecule check right and left side of the torsion and find largest
    * overlapping fragment. This can be used to align the optimized molecules.
    * @param coreFile name of output file
    */
   public void computeCore(String coreFile) throws IllegalArgumentException
   {  List<OEGraphMol> firstFrag =  new ArrayList<OEGraphMol>(cachedInputMols.size());
      List<OEGraphMol> secndFrag =  new ArrayList<OEGraphMol>(cachedInputMols.size());

      // loop over molecules, delete torsion bond,
      // store left side in firstFrag and right side in secndFrag
      for(OEMolBase mol: cachedInputMols)
      {  OEAtomBase[] molBondAts = getTorsionAtoms(mol, bondAtoms);
         OEBondBase bd = mol.GetBond(molBondAts[1],molBondAts[2]);
         int bondType = bd.GetOrder();
         mol.DeleteBond(bd);

         molBondAts[1].SetIntData(FRAGNumTag, 1); // remember atoms that where connecting to torsion
         molBondAts[2].SetIntData(FRAGNumTag, 2);

         int[] parts = new int[mol.GetMaxAtomIdx()];
         int pcount = oechem.OEDetermineComponents(mol, parts);
         OEPartPredAtom pred = new OEPartPredAtom(parts);
         for (int i = 1; i <= pcount; ++i)
         {  pred.SelectPart(i);
            OEGraphMol partmol = new OEGraphMol();
            oechem.OESubsetMol(partmol, mol, pred, true);

            OEAtomBaseIter atIt = partmol.GetAtoms();
            while(atIt.hasNext())
            {  OEAtomBase at1 = atIt.next();
               int fragNum = at1.GetIntData(FRAGNumTag);
               if( fragNum == 1 )
               {  OEGraphMol newFrag = new OEGraphMol(partmol);
                  at1 = AtomIntPropertyFunctor.getFirst(newFrag, FRAGNumTag, fragNum);
                  OEAtomBase at2 = newFrag.NewAtom(molBondAts[2]);
                  newFrag.NewBond(at1, at2, bondType);
                  firstFrag.add(new OEGraphMol(newFrag));
                  break;
               }
               if( fragNum == 2 )
               {  OEGraphMol newFrag = new OEGraphMol(partmol);
                  at1 = AtomIntPropertyFunctor.getFirst(newFrag, FRAGNumTag, fragNum);
                  OEAtomBase at2 = newFrag.NewAtom(molBondAts[1]);
                  newFrag.NewBond(at1, at2, bondType);
                  secndFrag.add(newFrag);
                  break;
               }
            }
            atIt.delete();
            partmol.delete();
         }
         pred.delete();
      }

      // find MCSS of all fragments
      CoreMatchResult coreMat = getMCSSCore(firstFrag);
      CoreMatchResult secndCoreMat = getMCSSCore(secndFrag);

      // find the one with more matches or more atoms
      if( coreMat.matchCount < secndCoreMat.matchCount
               || (coreMat.matchCount == secndCoreMat.matchCount
                   && coreMat.core.NumAtoms() < secndCoreMat.core.NumAtoms()))
      {  coreMat.close();
         coreMat = secndCoreMat;
      }else
      {  secndCoreMat.close();
      }

      for( OEMolBase mol : firstFrag) mol.delete();
      for( OEMolBase mol : secndFrag) mol.delete();

      OEGraphMol core = coreMat.getUnchargedCore();

      oemolothread ofs = new oemolothread(coreFile);
      oechem.OEWriteMolecule(ofs, core);
      ofs.close();
      ofs.delete();
      core.delete();
      coreMat.close();
   }

   /*
    * Run the torsion scanner which generates a set of conformers.  Assumes parameters have been set up through constructor
    * @inMol the input molecule\
    * @return a multi-conformer molecule
    */
   public OEMCMolBase run (OEMolBase inMol) throws IOException
   {
      OEMCMolBase torsionConformers = null;

      // Cache mol to use for later core file construction
      if (coreFilename != null)
         cachedInputMols.add(new OEGraphMol(inMol));

      try
      {
         // Get the atoms for this molecule based on coords of the
         //   bondFile atoms
         OEAtomBase[] torsionAtoms = getTorsionAtoms(inMol, bondAtoms);

         // Mark these atoms as "fixed"
         OEBondBase torBond = inMol.GetBond(torsionAtoms[1], torsionAtoms[2]);
         torBond.SetBoolData(ConformerSampler.ISFixedBondTag, true);

         // Set the fixed atom SDF tag
         setTorsionAtomSDTag(inMol, getTorsionAtomIndices(inMol, torsionAtoms));

         //
         // Main angular scan
         // Generate the fixed angle conformers around the torsion. (e.g. 0,10,20,30...)
         //
         torsionConformers = generateScannedTorsionConformers (inMol, torsionAtoms);

         // Possibly minimize or expand conformers
         if (this.nMaxConfsPerStep > 1)
         {
            // If nMaxConfsPerStep > 1,  generate more conformers at each rotatable bond outside the torsionAtoms bond.
            // These new confs will be written to torsionConformers
            expandConformers (torsionConformers, torsionAtoms, doMinimize);
         }
         else
         {
            if (doMinimize)
            {
               // Minimize each conformer
               minimizeConformers (torsionConformers, torsionAtoms);
            }
         }

         // Cleanup
         torsionAtoms = null;
         torBond = null;
      }
      catch(Exception e)
      {
         System.err.println(e.getMessage() + " compound ignored!");
      }

      return torsionConformers;
   }


   /*
    * Internal method to initialize the conformer sampler
    */
   private void initializeConformerSampler()
   {
      if (this.conformerSampler == null)
         this.conformerSampler = new ConformerSampler(TORSION_SCANNER_ANGLE_FILENAME);
   }

   /**
    * Internal method to read a file containing four atoms defining
    * the torsion bond.
    * @param bondFilename filename of sdf file or mol with four torsion atoms
    * @return a mol object containing the torsion atoms
    * @throws Error
    */
   private static OEMolBase readBondAtomsFromFile(String bondFilename) throws Error
   {
      oemolithread ifs = new oemolithread(bondFilename);
      OEGraphMol   mol = new OEGraphMol();

      try
      {
         oechem.OEReadMolecule(ifs, mol);
         if (oechem.OEReadMolecule(ifs, mol))
            throw new Error("BondMol has more than one record");
      }
      finally
      {
         ifs.close();
         ifs.delete();
      }

      if (!mol.IsValid())
         throw new Error("BondMol is inValid");

      if (mol.NumAtoms() != 4)
         throw new Error("BondMol does not have 4 atoms");

      return mol;
   }

   /**
    * Initialize SZYBKI as a minimizer
    * @param doMinimize true if a minimization should be used
    * @param constraintStrength
    * @return a szybki instance
    */
   private OESzybki initializeSZYBKI(String constraintStrength)
   {
      OESzybki szybki = null;

      szybki = new OESzybki(OEForceFieldType.MMFF94S);
      szybki.SetSolventModel(OESolventModel.Sheffield);
      // discussed with Ignacio, lets do No minimization, full minimization
      // constrained minimization
      // with default on full minimization

      double forceconstant = 90; // 90 is strong leaving only little space
                                 // for movement
      double distance = 0;
      if ("strong".equalsIgnoreCase(constraintStrength)
               || constraintStrength == null
               || constraintStrength.length() == 0)
         forceconstant = 90;
      else if ("medium".equalsIgnoreCase(constraintStrength))
         forceconstant = 45;
      else if ("weak".equalsIgnoreCase(constraintStrength))
         forceconstant = 20;
      else if (!"none".equalsIgnoreCase(constraintStrength))
         forceconstant = Double.parseDouble(constraintStrength);

      if (!"none".equalsIgnoreCase(constraintStrength))
         szybki.SetHarmonicConstraints(forceconstant, distance,
                  new OEIsHeavy());

      return szybki;
   }


   /*
    * Given a mol that is presumed to contain 4 atoms, reorder the atoms
    * @param mol a molecule with four bonded torsion atoms
    * @return ordered set of atoms
    */
   private static OEAtomBase[] getOrderedBondAtoms(OEMolBase mol)
   {
      OEAtomBaseIter atIt = mol.GetAtoms();
      OEAtomBase [] atoms = new OEAtomBase[mol.NumAtoms()];
      OEAtomBase startAtom = null;
      int nAt = 0;
      while (atIt.hasNext())
      {  OEAtomBase at = atIt.next();
         int deg = at.GetExplicitDegree();
         if (deg == 1 && startAtom == null)
         {  startAtom = at;
         }

         if (deg > 2)
            throw new Error(String.format(
                     "Atom %d has more more (%d) than two bonds", nAt, deg));

         atoms[nAt++] = at;
      }
      atIt.delete();
      if (startAtom == null)
         throw new Error("No atom whit one bond found!");

      // use the single bonded atom to traverse to the four connected atoms and
      // define the atoms in sequence
      atoms[0] = startAtom;
      atoms[1] = getOtherAtom(atoms[0]);
      atoms[2] = getOtherAtom(atoms[1], atoms[0]);
      atoms[3] = getOtherAtom(atoms[2], atoms[1]);

      return atoms;
   }





   /*
    * Minimize a set of conformers
    * @param torsionConformers a set of torsion conformers to be minimized
    * @param torsionAtoms atoms to be held fixed
    */
   private void minimizeConformers(OEMCMolBase torsionConformers, OEAtomBase[] torsionAtoms)
   {
      System.err.println(String.format("%d Conformers will be minimized at each step.  \n\tMinimizing...", torsionConformers.NumConfs()));

      for (OEConfBase torsionConf : torsionConformers.GetConfs())
      {
         oechem.OECopySDData(torsionConf, torsionConformers);
         minimize(torsionConf, torsionAtoms);
      }
   }


   /*
    * Expand a set of conformers.  the input set will be changed to contain the expanded set
    * @torsionConformers an input set of conformers to be expanded
    * @torsionAtoms  atoms to be held fixed
    *
    */
   private void expandConformers(OEMCMolBase torsionConformers, OEAtomBase[] torsionAtoms, boolean doMinimize)
   {
      OEMCMolBase generatedConfs = null;

      // Initialize conformer sampler if needed
      initializeConformerSampler();

      if (doMinimize) {System.err.println(String.format("%d additional conformers per step will be minimized to a single conformer per step. \n\tMinimizing...", nMaxConfsPerStep));}

      // Copy the metadata
      OEMCMolBase additionalConfs = new OEMol(torsionConformers);
      additionalConfs.DeleteConfs();

      OEConfBaseIter confIt = torsionConformers.GetConfs();
      while( confIt.hasNext() )
      {  OEConfBase torsionConf = confIt.next();
         // Make sure the OEConfBase contains the SD data
         oechem.OECopySDData(torsionConf, torsionConformers);

         generatedConfs = conformerSampler.createConformations(torsionConf, nMaxConfsPerStep);

         if (doMinimize)
         {
            // Minimize the original plus expanded conformers.  Only keep the lowest energy single conformer

            // Add the original
            generatedConfs.NewConf(torsionConf);

            OEConfBase singleLowEMol = getLowestEMinimizedConf (generatedConfs, torsionAtoms);
            addConformer(additionalConfs, singleLowEMol);
         }
         else
         {
            // Add original conf
            additionalConfs.NewConf(torsionConf);

            // Add all unminimized generated confs
            addMultiConformers(additionalConfs, generatedConfs);
         }
      }
      confIt.delete();

      // Finally add the stored confs to the passed in set.  Since we read saved the
      // originals, delete what was there and add them back to preserve order.
      torsionConformers.DeleteConfs();
      addMultiConformers(torsionConformers, additionalConfs);

   }




   /*
    * Initialize the internal multi-conf molecule using the passed in molecule
    * @param inMol an input molecule to use for initialization
    * @return a multi-conformer molecule
    */
   private static OEMCMolBase initMCMol(OEMolBase inMol)
   {
      OEMCMolBase torsionConformers = null;
      if (inMol == null)
      {
         throw new Error("Internal error, no input molecule.");
      }
      else
      {
         torsionConformers = new OEMol(inMol);
      }
      return torsionConformers;
   }


   /*
    * Add a conformer to the destination OEMCMol
    * @param destMCMol the destination MCMol
    * @param conformer the input conformer
    */
   private static void addConformer(OEMCMolBase destMCMol, OEConfBase conformer)
   {
      if (destMCMol == null)
      {
         throw new Error("Internal error, trying to add conformer to null object");
      }

      destMCMol.NewConf(conformer);
   }


   /*
    * Add MCMolBase to a destination destMCMol
    * @param destMCMol the destination multi conformer molecule to which confs will be added
    * @param confsToAdd the conformers to add
    */
   private static void addMultiConformers(OEMCMolBase destMCMol, OEMCMolBase confsToAdd)
   {
      for (OEConfBase conf : confsToAdd.GetConfs())
      {
         addConformer(destMCMol, conf);
      }
   }


   /*
    * Get the atom indices for the torsion atoms in the incoming molecule
    * @param mol the incomin molecule
    * @param torsionAtoms  array of atoms whose indices are to be returned
    * @return an array of indices
    */
   private static int [] getTorsionAtomIndices(OEMolBase mol, OEAtomBase[] torsionAtoms)
   {
      int[] torsionAtomIndices = new int[4];
      int nTorsionIndex = 0;

      // Searching to maintain order of torsion Atoms
      for (OEAtomBase torsionAtom : torsionAtoms)
      {
         int nMolAtomIndex = 0;
         for (OEAtomBase molAtom : mol.GetAtoms())
         {

            if (molAtom.GetIdx() == torsionAtom.GetIdx())
            {
               torsionAtomIndices[nTorsionIndex]=nMolAtomIndex;
               nTorsionIndex++;
            }
            nMolAtomIndex++;
         }
      }
      return torsionAtomIndices;
   }


  /*
   * Write out the atom indices which define the torsion
   * These should match the index of the atom in the written out sdf file,
   *  says OpenEye
   *  @param mol the input molecule
   *  @param torsionAtomIndices the indices to write out
   */
   private void setTorsionAtomSDTag(OEMolBase mol, int[] torsionAtomIndices)
   {
      String torsionAtomsTagValue = String.format("%d %d %d %d",
               torsionAtomIndices[0], torsionAtomIndices[1], torsionAtomIndices[2], torsionAtomIndices[3]);

      oechem.OESetSDData(mol, this.torsionAtomsTag, torsionAtomsTagValue);
   }

   /*
    * Get the single lowest energy conformer from a set of input conformers.  Minimize each one
    * @param mcMol the input multi-conformer mol which will be reduced to a single conformer
    * @param fixAtoms an array of atoms to be held fixed during minimization
    * @return a single conformer molecule
    */
   private OEConfBase getLowestEMinimizedConf(OEMCMolBase mcMol, OEAtomBase[] fixAtoms)
   {
      double minE = Double.MAX_VALUE;
      OEConfBase lastMinEConf  = null;
      OEMCMolBase minConfMol = new OEMol(mcMol);

      minConfMol.DeleteConfs();

      OEConfBaseIter confIt = mcMol.GetConfs();
      while( confIt.hasNext() )
      {  OEConfBase currConf = confIt.next();

         double newE = minimize(currConf, fixAtoms);
         if (minE > newE)
         {
            // We have a new low energy
            minE = newE;
            lastMinEConf = currConf;
         }
      }
      confIt.delete();

      return lastMinEConf;
   }


   /**
    * Internal method to minimize a conformer and return the energy
    * @param newConf
    * @param fixAtoms atoms to hold fixed
    * @return energy of the minimized conformer
    */
   private double minimize (OEConfBase newConf, OEAtomBase[] fixAtoms)
   {  AtomListFunctor atFunct = new AtomListFunctor(fixAtoms);

      if (!szybki.FixAtoms(atFunct))
      {  atFunct.delete();
         throw new Error("Failed to fix atoms");
      }
      OEMol copy = new OEMol(newConf);
      oechem.OECopySDData(copy, newConf);

      double res = Double.MAX_VALUE;
      OESzybkiResults szRes = new OESzybkiResults();
      if (szybki.call(newConf, szRes))
          res = szRes.GetFinalTotalPotential();

      szRes.delete();
      oechem.OECopySDData(newConf, copy);
      copy.delete();
      szybki.ClearFixAtoms();
      atFunct.delete();

      return res;
   }

   /*
    * Find first bonded atom to an input atom
    * @param at input atom
    * @return the "next" bonded atom
    */
   private static OEAtomBase getOtherAtom(OEAtomBase at)
   {
      OEAtomBaseIter atIt = at.GetAtoms();
      OEAtomBase at2 = null;
      if (atIt.hasNext())
         at2 = atIt.next();

      assert !atIt.hasNext() : "Atom has more than one bond: " + at.GetIdx();
      atIt.delete();
      return at2;
   }

   /*
    * Get first atom attached to at that is not excludedAt
    * @param at the input atom
    * @param excludedAt an atom that we don't want to return
    * @return an atom bonded to the input and not  excludedAt
    */
   private static OEAtomBase getOtherAtom(OEAtomBase at, OEAtomBase excludedAt)
   {
      int exclIdx = excludedAt.GetIdx();
      OEAtomBaseIter atIt = at.GetAtoms();
      OEAtomBase at2 = null;
      while (atIt.hasNext())
      {
         OEAtomBase oa = atIt.next();
         if (oa.GetIdx() != exclIdx)
         {
            at2 = oa;
            break;
         }
      }
      atIt.delete();

      return at2;
   }



   /** find MCSS of molecules **/
   private static CoreMatchResult getMCSSCore(List<OEGraphMol> molList)
   {  if( molList.size() == 0 ) return null;

      OEGraphMol core = new OEGraphMol(molList.get(molList.size()-1));
      oechem.OESuppressHydrogens(core);
      if( molList.size() == 1 ) return new CoreMatchResult(core,1);

      int atomexpr = 0; // ignore atom type
      int bondexpr = 0; // ignore bond order
      int mcsstype = OEMCSType.Exhaustive;
      OEMCSSearch mcss = new OEMCSSearch(core, atomexpr, bondexpr, mcsstype);
      mcss.SetMCSFunc(new OEMCSMaxBondsCompleteCycles());
      mcss.SetMinAtoms(3);
      mcss.SetMaxMatches(1);
      int matchCount=0;
      // start from the back so that we are left we coordinates of first
      for(int i=molList.size()-2; i>=0; i--)
      {  OEGraphMol mol = molList.get(i);
         oechem.OESuppressHydrogens(mol);
         OEMatchBaseIter mIt = mcss.Match(mol);
         if(! mIt.hasNext())
         {  mIt.delete();
            continue;
         }
         matchCount++;
         core.Clear();
         oechem.OESubsetMol(core,mIt.next(),true);

         // reinitialize to smaller core
         mcss.Init(core, atomexpr, bondexpr);
         mIt.delete();
      }
      mcss.delete();

      return new CoreMatchResult(core,matchCount);
   }

   /**
    * Compute nSteps conformations starting from mol rotating around the torsion
    * defined by molBondAts
    * @param startMol input molecule
    * @param molBondAts 4 atoms defining the torsion around the two central atoms
    * @return a multiconf mol
    */
   private OEMCMolBase generateScannedTorsionConformers (OEMolBase startMol, OEAtomBase[] molBondAts)
   {
      // Creating new
      OEMCMolBase confMol = initMCMol(startMol);

      // if startTorsion was not supplied (NaN) use torsion in input mol
      double myStartTorsion = nStartTor;

      if( Double.isNaN(nStartTor) )
         myStartTorsion = Math.toDegrees(oechem.OEGetTorsion(
                  confMol,molBondAts[0],molBondAts[1], molBondAts[2],molBondAts[3]));

      OEConfBaseIter confIt = confMol.GetConfs();
      OEConfBase baseConf = confIt.next();
      confIt.delete();

      for (int s = 0; s < nSteps; s++)
      {
         OEConfBase newConf = confMol.NewConf(baseConf);
         oechem.OEAddExplicitHydrogens(newConf, false, true);
         double angle = myStartTorsion + s * this.nTorIncr;

         oechem.OESetSDData( newConf, ANGLE_TAG, String.format("%.1f", angle));

         oechem.OESetTorsion(newConf, molBondAts[0],molBondAts[1],
                                      molBondAts[2],molBondAts[3], Math.toRadians(angle));

      }

      // Now delete the starting molecule which was passed in
      confMol.DeleteConf(baseConf);

      return confMol;
   }


   /**
    * get atoms closest to bondAtoms in mol
    */
   private OEAtomBase[] getTorsionAtoms(OEMolBase mol, OEAtomBase[] bondAtoms) throws IllegalArgumentException
   {  OEAtomBase[] molBondAts = new OEAtomBase[bondAtoms.length];

      for (int i = 0; i < bondAtoms.length; i++)
      {
         double minDist = Double.MAX_VALUE;
         OEAtomBaseIter atIt = mol.GetAtoms();

         while (atIt.hasNext())
         {
            OEAtomBase at = atIt.next();
            double dist = oechem.OEGetDistance(mol, at, bondMol, bondAtoms[i]);
            if (dist < minDist)
            {
               molBondAts[i] = at;
               minDist = dist;
            }
         }
         atIt.delete();
      }

      for (int i = 0; i < molBondAts.length; i++)
         for (int j = i+1; j < molBondAts.length; j++)
            if (molBondAts[i].GetIdx() == molBondAts[j].GetIdx())
               throw new IllegalArgumentException(
                        String.format(
                                 "Atom %d is nearest neighbor to two atoms in bond file %d and %d for %s",
                                 molBondAts[i].GetIdx(), bondAtoms[i].GetIdx(),
                                 bondAtoms[j].GetIdx(), OETools.molToCanSmi(mol, true)));

      OEBondBase bd = mol.GetBond(molBondAts[1], molBondAts[2]);
      if( bd == null ) throw new Error("No bond between second and third atom");
      if( bd.IsInRing() ) throw new Error("Central bond is in Ring");
      if( bd.GetIntType() != 1 ) throw new Error("Central bond must be single bond");

      return molBondAts;
   }
}

/*
 * Class containing a Max Common Substructure search with the number of matches
 */
class CoreMatchResult
{  OEGraphMol core;
   int matchCount;

   /**
    * Constructor
    * @param core molecule object match
    * @param matchCount number of matches
    */
   public CoreMatchResult(OEGraphMol core, int matchCount)
   {  oechem.OESuppressHydrogens(core);
      this.core = core;
      this.matchCount = matchCount;
   }

   public OEGraphMol getUnchargedCore()
   {  OEGraphMol mcopy = new OEGraphMol(core);
      OEAtomBaseIter atit = mcopy.GetAtoms();
      while(atit.hasNext())
      {  OEAtomBase at = atit.next();
         at.SetFormalCharge(0);
      }
      atit.delete();
      return mcopy;
   }
   public void close() { core.delete(); }
}
