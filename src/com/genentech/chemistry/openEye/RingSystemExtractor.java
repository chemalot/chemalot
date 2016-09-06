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
package com.genentech.chemistry.openEye;


import openeye.oechem.OEAtomBase;
import openeye.oechem.OEBondBase;
import openeye.oechem.OEGraphMol;
import openeye.oechem.OEMolBase;
import openeye.oechem.oechem;

import com.aestel.utility.Message;
import com.aestel.utility.exception.IncorrectInputException;
import com.genentech.oechem.tools.OETools;


/**
 * Use it to extract the largest linked ring system and the basic ring systems 
 * of a molecule.
 *
 * @author Man-Ling Lee / July 31, 2010
 * Copyright 2010 Genentech
 */
public final class RingSystemExtractor
{  private OEMolBase compoundMol;
   private OEMolBase largestRingSystem;
   private OEMolBase basicRingSystems;

   
   public RingSystemExtractor()
   {  compoundMol       = new OEGraphMol();
      largestRingSystem = new OEGraphMol();
      basicRingSystems  = new OEGraphMol();
   }
   
   
   public boolean hasLargestRingSystem()
   {  if( largestRingSystem.GetAtoms().hasNext() )
         return true;
      return false;
   }
   
   public String getLargestRingSystemSMILES()
   {  return OETools.molToCanSmi( largestRingSystem, true); }
   
   
   public int getBasicRingSystemCount()
   {  String smi = OETools.molToCanSmi( basicRingSystems, true );
      if( smi.length() == 0 )
         return 0;
      String[] fs = smi.split( "\\." );
      return fs.length;
   }
   
   public String getBasicRingSystemsSMILES()
   {  return OETools.molToCanSmi( basicRingSystems, true ); }
   
   
   /**
    * Extract the largest ring system and the basic ring systems from smiles.
    */
   public void extractFromSMILES( String smiles )
   throws IncorrectInputException
   {  compoundMol.Clear();
      if( oechem.OEParseSmiles( compoundMol, smiles ) )
      {  extract();
      } else
      {  String s = "SMILES string was invalid!" ;
         Message msg = new Message( s, Message.Level.ERROR, null );
         throw new IncorrectInputException( msg );
      }
   }
   
   /**
    * Extract the largest ring system and the basic ring systems from mol.
    * 
    * Note that a copy of mol will be made.
    */
   public void extract( OEMolBase mol )
   {  compoundMol.Clear();
      oechem.OEAddMols( compoundMol, mol );
      extract();
   }

   private void extract()
   {  //TODO: Check if the stereochemistry of a heavy atom is preserved if 
      //      the stereochemistry is express by the bond to the explicit 
      //      hydrogen. Currently retainStereo (#3) is set to false.
      oechem.OESuppressHydrogens(compoundMol, false, false, false );
      oechem.OEAssignAromaticFlags( compoundMol );
      
      extractLargestRingSystem();
      extractbasicRingSystems();
      
   }
   
   /**
    * Delete terminal chains from the molecule.
    * 
    * Exo-cyclic double bond are considered as part of the ring.
    */
   private void extractLargestRingSystem()
   {  largestRingSystem.Clear();
      oechem.OEAddMols( largestRingSystem, compoundMol );
      for( OEAtomBase atom : largestRingSystem.GetAtoms() )
        deleteSideChainAtoms( atom );
   }
   
   
   /**
    * Delete the atoms in a side chain.
    * 
    * If atom is a terminal atom, all the atoms in this side chain will be
    * deleted from largestRingSystem.
    */
   private void deleteSideChainAtoms( OEAtomBase atom )
   {  if( atom.IsInRing() )
         return;
      
      //Get the neighbor atoms of the current atom and.
      OEAtomBase neigborAtom = null;
      OEBondBase neigborBond = null;
      int neigborAtomCount = 0;
      for( OEBondBase bond : atom.GetBonds() )
      {  ++neigborAtomCount;
         if( neigborAtomCount == 1 )
         {  neigborBond = bond;
            neigborAtom = bond.GetNbr( atom );
         } else
            break;
      }
      if( neigborAtomCount == 0 )    //Single atom, e.g. O for water
      {  largestRingSystem.DeleteAtom( atom );
         return;
      }
      
      //Stop here if atom is not a non-terminal atom and an atom at one end of 
      //the exo-cyclic double bond
      if( neigborAtomCount > 1 ||
          ( neigborBond.GetOrder() > 1 && neigborAtom.IsInRing() ) )
         return;
      
      largestRingSystem.DeleteAtom( atom );
      neigborAtom.SetImplicitHCount( 
                  neigborAtom.GetImplicitHCount() + neigborBond.GetOrder() );
      deleteSideChainAtoms( neigborAtom );
   }

   
   /**
    * Delete chains linking the basic ring system.
    * 
    * Exo-cyclic double bond are considered as part of the ring.
    */
   private void extractbasicRingSystems()
   {  basicRingSystems.Clear();
      oechem.OEAddMols( basicRingSystems, largestRingSystem );
      
      //Remove linker bonds
      for( OEBondBase bond : basicRingSystems.GetBonds() )
      {  if( bond.IsInRing() )
            continue;
      
         //Keep the exo-cyclic double bond.
         OEAtomBase beginAtom = bond.GetBgn();
         OEAtomBase endAtom   = bond.GetEnd();
         if( bond.GetOrder() > 1 )
         {  if( beginAtom.IsInRing() || endAtom.IsInRing() )
               continue;
         }
         beginAtom.SetImplicitHCount( 
                  beginAtom.GetImplicitHCount() + bond.GetOrder() );
         endAtom.SetImplicitHCount( 
                  endAtom.GetImplicitHCount() + bond.GetOrder() );
         basicRingSystems.DeleteBond( bond );
      }
      //Remove the linker atoms
      for( OEAtomBase atom : basicRingSystems.GetAtoms() )
      {  if( !atom.GetBonds().hasNext() )
            basicRingSystems.DeleteAtom( atom );
      }
   }
   
   
   public static void main(String [] args) 
   throws IncorrectInputException
   {
      String smi1 = "C1CCCC1=C(CC(N)C)Cc2cc(ccc2)C(=NC)C";
      RingSystemExtractor rsExtractor = new RingSystemExtractor();
      
      rsExtractor.extractFromSMILES( smi1 );
      
      System.out.println( "Input:\t" + smi1 + "\n\n" );
      System.out.println( "largest RS\t" + rsExtractor.getLargestRingSystemSMILES() + "\n\n" );
      System.out.println( "basic RSs\t" + rsExtractor.getBasicRingSystemsSMILES() + "\n\n" );
      
      
      
   }
}
