/*
   Copyright 2008-2014 Genentech Inc.

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
package com.genentech.oechem.tools;

import openeye.oechem.*;

public class OETools
{
   private static final AtomHasChargeFunctor HASChargeFunctor = new AtomHasChargeFunctor();

   /** get the molfile of this OEMolBase */
   public static String molToString(OEMolBase mol)
   {  return molToString(mol, OEFormat.MDL);
   }

   public static String molToSDF(OEMolBase mol)
   {  return molToString(mol, OEFormat.SDF);
   }

   /** get the molfile of this OEMolBase */
   private static String molToString(OEMolBase mol, int type)
   {
      oemolostream ofs = new oemolostream();
      ofs.openstring();
      ofs.SetFormat(type);
      OEGraphMol myMol = new OEGraphMol(mol); // OEWriteMolecule will modify
                                              // mol!
      oechem.OEWriteMolecule(ofs, myMol);
      String molStr = ofs.GetString();
      ofs.close();
      ofs.delete();
      myMol.delete();

      return molStr;
   }

   /**
    * Convert a string containing an isomeric smiles to a mol object. The mol
    * object is not cleared before loading the new molecule.
    *
    * @return the same mol object as passed in on input, just for convenience
    */
   public static OEMolBase smiToMol(OEMolBase mol, String smi)
   {
      if (smi == null || smi.length() == 0)
         return mol;

      oemolistream ifs = new oemolistream();
      ifs.openstring(smi);
      ifs.SetFormat(OEFormat.ISM);

      oechem.OEReadMolecule(ifs, mol);
      ifs.close();
      ifs.delete();

      return mol;
   }

   /**
    * convert a string containing an MDL molfile to a mol object.
    *
    * @return the same mol object as passed in on input, just for convenience
    */
   public static OEMolBase stringToMol(OEMolBase mol, String molStr)
   {
      oemolistream ifs = new oemolistream();
      ifs.openstring(molStr);
      ifs.SetFormat(OEFormat.MDL);

      oechem.OEReadMolecule(ifs, mol);
      ifs.close();
      ifs.delete();

      return mol;
   }

   /**
    * This method might change internal structure of the molecule by aromatizing
    * rings etc.
    *
    * @param mol
    *           molecule to convert to smiles
    * @param iso
    *           if true the isomeric smiles is created
    *
    * @return Canonical smiles
    */
   public static String molToCanSmi(OEMolBase mol, boolean iso)
   {

      // implementation could use oechem.OECreateIsoSmiString(mol)
      // but beware that you have to consult the documentation about preparing
      // mol for that (eg. aromatization)

      int format = OEFormat.ISM;
      if (!iso)
         format = OEFormat.CAN;

      oechem.OESuppressHydrogens(mol, false, false, true);
      oechem.OEPerceiveChiral(mol);
      if( mol.GetDimension() > 2 && oechem.OEGetDimensionFromCoords(mol) <= 2)
         mol.SetDimension(2);


      oemolostream ofs = new oemolostream();
      ofs.openstring();
      ofs.SetFormat(format);
      oechem.OEWriteMolecule(ofs, mol);
      String molStr = ofs.GetString().trim();
      ofs.close();
      ofs.delete();

      // remove any trailing identifiers
      int sPos = molStr.indexOf(" ");
      if (sPos >= 0)
         molStr = molStr.substring(0, sPos);

      return molStr;
   }

   // public static final void SuppressHydrogens(OEMolBase oeMol) {
   // for ( OEAtomBaseIter aIt = in.GetAtoms(); aIt.hasNext(); ) {
   // OEAtomBase at = aIt.next();
   // if(! oechem.OEHasAtomStereoHydrogens(at))
   // at.;
   //
   // }

   /**
    * Find the average bond length in this molecule.
    */
   public static double getAverageBondLength(OEMolBase in)
   {
      double avrgBond = 0;
      OEBondBaseIter bIt = in.GetBonds();
      while (bIt.hasNext())
      {
         OEBondBase bd = bIt.next();
         OEAtomBase a1 = bd.GetBgn();
         OEAtomBase a2 = bd.GetEnd();

         avrgBond += oechem.OEGetDistance(in, a1, a2);
      }
      bIt.delete();
      avrgBond /= in.NumBonds();
      return avrgBond;
   }

   /**
    * Add toAdd to the base mol by looking at the 2D coordinates and orienting
    * toAdd to the right (bigger x) coordinate of base.
    *
    * @param base
    *           will be modified in place.
    * @param toAdd
    *           will be added to base will be translated in the process.
    * @return modified base.
    */
   public static OEMolBase combineStructures(OEMolBase base, OEMolBase toAdd)
   {  return combineStructures(base, toAdd, 1);
   }

   /**
    * @param count
    *           number of repetition of toAdd to be added
    * @param toAdd
    */
   public static OEMolBase combineStructures(OEMolBase base, OEMolBase toAdd,
            int count)
   {  assert base.IsValid() : "base molecule is invalid";
      assert toAdd.IsValid() : "added molecule is invalid";

      double baseMaxX = getMaxCoor(base, 0);
      double baseMinY = getMinCoor(base, 1);
      double addMinX = getMinCoor(toAdd, 0);
      double addMaxY = getMaxCoor(toAdd, 1);
      double addheight = addMaxY - getMinCoor(toAdd, 1);

      double avrgBondLength = OETools.getAverageBondLength(base);
      if (avrgBondLength <= 0.001 || Double.isNaN(avrgBondLength))
         avrgBondLength = 1;

      // X translation + over compensate for one Y shift
      oechem.OETranslate(toAdd, new double[]
         { baseMaxX + avrgBondLength / 1D - addMinX,
                  baseMinY - addMaxY - avrgBondLength / 1.8d, 0 });

      for (int i = 0; i < count; i++)
      { // Y translation
         oechem.OETranslate(toAdd, new double[]
            { 0, +addheight + avrgBondLength / 1.8d, 0 });

         oechem.OEAddMols(base, toAdd);
      }
      return base;
   }

   public static double getMinCoor(OEMolBase toAdd, int axis)
   {
      double min = Double.POSITIVE_INFINITY;
      double[] coor = new double[3];
      OEAtomBaseIter aIt = toAdd.GetAtoms();
      while (aIt.hasNext())
      {
         OEAtomBase at = aIt.next();
         toAdd.GetCoords(at, coor);
         if (coor[axis] < min)
            min = coor[axis];
      }
      aIt.delete();

      if (min == Double.POSITIVE_INFINITY)
         return 0; // no atoms
      return min;
   }

   public static double getMaxCoor(OEMolBase toAdd, int axis)
   {
      double max = Double.NEGATIVE_INFINITY;
      double[] coor = new double[3];
      OEAtomBaseIter aIt = toAdd.GetAtoms();
      while (aIt.hasNext())
      {
         OEAtomBase at = aIt.next();
         toAdd.GetCoords(at, coor);
         if (coor[axis] > max)
            max = coor[axis];
      }
      aIt.delete();

      if (max == Double.NEGATIVE_INFINITY)
         return 0; // no atoms
      return max;
   }

   /** get molecular weight using most abundant isotope of unspecified atoms */
   public static double getIsotpictWeight(OEMolBase mol)
   {
      double weight = 0D;
      OEAtomBaseIter aIt = mol.GetAtoms();
      while (aIt.hasNext())
      {
         OEAtomBase atom = aIt.next();

         int elemno = atom.GetAtomicNum();
         int mass = atom.GetIsotope();
         int implicitH = atom.GetImplicitHCount();
         if (elemno != 0 && mass != 0 && oechem.OEIsCommonIsotope(elemno, mass))
            weight += oechem.OEGetIsotopicWeight(elemno, mass);
         else
            weight += oechem.OEGetIsotopicWeight(elemno,oechem.OEGetDefaultMass(elemno));

         weight += implicitH * oechem.OEGetIsotopicWeight(OEElemNo.H, 1);
      }
      aIt.delete();

      return weight;
   }

   public static int getCharge(OEMolBase mol)
   {  OEAtomBaseIter atIt = mol.GetAtoms();
      int sum = 0;
      while( atIt.hasNext() )
      {  OEAtomBase at = atIt.next();
         sum +=  at.GetFormalCharge();
      }
      atIt.delete();

      return sum;
   }

   public static int countAtoms(OEMolBase mol, OEUnaryAtomPred pred)
   {  int counts = 0;

      OEAtomBaseIter atIt = mol.GetAtoms();
      while( atIt.hasNext() )
         if (pred.constCall(atIt.next())) counts++;
      atIt.delete();

      return counts;
   }

   public static void neutralize(OEMolBase mol)
   {  OEAtomBaseIter atIt = mol.GetAtoms(HASChargeFunctor);
      while(atIt.hasNext())
      {  OEAtomBase at = atIt.next();
         if( at.GetFormalCharge() < 0 )
         {  neutralizeNegative(at);

         } else if( at.GetFormalCharge() > 0 )
         {  neutralizePositive(at);
         }
      }
      atIt.delete();
   }

   private static void neutralizeNegative(OEAtomBase at)
   {  assert at.GetFormalCharge() < 0;

      if(! (at.IsNitrogen() || at.IsOxygen() || at.IsSulfur() || at.IsPhosphorus()
          ||at.IsHalogen() ) )
         return;

      boolean hasChargedNeighbor = false;
      OEAtomBaseIter nIt = at.GetAtoms();
      while( nIt.hasNext() )
      {  if( nIt.next().GetFormalCharge() > 0 )
         {  hasChargedNeighbor = true;
            break;
         }
      }
      nIt.delete();

      if( hasChargedNeighbor ) return;

      at.SetImplicitHCount(at.GetImplicitHCount()+1);
      at.SetFormalCharge(at.GetFormalCharge()+1);
   }

   private static void neutralizePositive(OEAtomBase at)
   {  assert at.GetFormalCharge() > 0;

      if( at.GetTotalHCount() == 0 ||
         ! ( at.IsNitrogen() || at.IsPhosphorus() ) )
         return;

      boolean hasChargedNeighbor = false;
      OEAtomBase hAt = null;

      OEAtomBaseIter nIt = at.GetAtoms();
      while( nIt.hasNext() )
      {  OEAtomBase nAt = nIt.next();
         if( nAt.GetFormalCharge() < 0 )
         {  hasChargedNeighbor = true;
            break;
         }
         if( nAt.IsHydrogen() )
            hAt = nAt;
      }
      nIt.delete();

      if( hasChargedNeighbor ) return;

      if( at.GetImplicitHCount() > 0 )
      {  at.SetImplicitHCount(at.GetImplicitHCount()-1);
      } else
      {  assert hAt != null;
         at.GetParent().DeleteAtom(hAt);
      }
      at.SetFormalCharge(at.GetFormalCharge()-1);
   }


   public static void main(String... args)
   {  OEGraphMol mol = new OEGraphMol();
      OETools.smiToMol(mol, "[K]");
      System.err.println(OETools.getIsotpictWeight(mol));
      mol.delete();
   }

   public static void main2(String... args)
   {
      oemolistream ifs = new oemolistream(args[0]);
      oemolostream ofs = new oemolostream(args[1]);

      OEGraphMol mol1 = new OEGraphMol();
      OEGraphMol mol1Copy = new OEGraphMol();
      OEGraphMol mol2 = new OEGraphMol();
      oechem.OEReadMolecule(ifs, mol1);

      while (oechem.OEReadMolecule(ifs, mol2))
      {
         mol1Copy.Clear();
         oechem.OEAddMols(mol1Copy, mol1);
         oechem.OEWriteMolecule(ofs, combineStructures(mol2, mol1Copy, 5));
         mol2.Clear();
      }
      mol2.delete();
      mol1Copy.delete();
      mol1.delete();
      ofs.close();
      ofs.delete();
      ifs.close();
      ifs.delete();
   }
}
