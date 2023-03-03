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
package com.genentech.chemistry.openEye.aLogP;

import java.io.*;
import java.util.*;

import openeye.oechem.*;

import com.genentech.oechem.tools.Atom;




/**
 * Calculate ALogP value and counts of ALogP atom types as per:
 *
 * Ghose, A. K.; Viswanadhan, V. N.; Wendoloski, J. J.
 * Prediction of Hydrophobic (Lipophilic) Properties of Small Organic Molecules Using Fragmental Methods:
 * An Analysis of ALOGP and CLOGP Methods.
 * J. Phys. Chem. A 1998, 102, 3762-3772.

 * @author Alberto Gobbi/ 2012
 * Copyright 2012 Genentech
 */
public class ALogPCalculator
{  public static final int ALOGPAtIdxTag = oechem.OEGetTag("ALOGPAtIdxTag");

   private final Map<String,List<ALogPAtom>> atomTypesByAtomSym
              = new HashMap<String,List<ALogPAtom>>(10);
   private final List<ALogPAtom> allAtomTypes = new ArrayList<ALogPAtom>(130);
   private final int[] atomTypeCounts;
   private final OEGraphMol noExpHMol;
   private final OEGraphMol allExpHMol;
   private final boolean validatedAssignment;

   private boolean[] atomWasAssigned;
   private boolean[] hWasAssigned;


   public ALogPCalculator( String smartFile, boolean validatedAssignment )
   {  int maxAtomType = readSMARTS(smartFile);
      atomTypeCounts = new int[maxAtomType+1];

      this.validatedAssignment = validatedAssignment;
      noExpHMol = new OEGraphMol();
      allExpHMol = new OEGraphMol();
   }


   public void close()
   {  for(ALogPAtom  a : allAtomTypes)
         a.close();

      allExpHMol.delete();
      noExpHMol.delete();
   }

   public double computeALogP(OEMolBase mol)
   {  Arrays.fill(atomTypeCounts, 0);
      noExpHMol.Clear();
      allExpHMol.Clear();

      if( validatedAssignment )
      {  atomWasAssigned = new boolean[mol.GetMaxAtomIdx()+1];
         hWasAssigned    = new boolean[mol.GetMaxAtomIdx()+1];
         assignIndexTag(mol);
      }

      oechem.OEAddMols(allExpHMol, mol, (String)null);
      oechem.OEAddMols(noExpHMol, mol, (String)null);
      oechem.OEAddExplicitHydrogens(allExpHMol);
      oechem.OESuppressHydrogens(noExpHMol, false, false,false);

      double aLogP = 0D;
      aLogP += processAtoms(noExpHMol, false);
      aLogP += processAtoms(allExpHMol, true);

      if( validatedAssignment )
      {  validateAssignments(mol);
      }
      return aLogP;
   }


   private double processAtoms(OEMolBase mol, boolean hAreExplicit)
   {  double aLogP = 0D;

      OEAtomBaseIter atIt = mol.GetAtoms();
      while( atIt.hasNext() )
      {  OEAtomBase at = atIt.next();

         if( at.GetAtomicNum() == 1 ) continue; // H are considered on central atom

         aLogP += getAtomHContrib( at, hAreExplicit );
         //System.err.printf("Hcontrib %f\n",getAtomHContrib( at, hAreExplicit ));

         String atSym = oechem.OEGetAtomicSymbol(at.GetAtomicNum());

         if( "H".equals(atSym) ) continue;  // H's are considered on parent atoms

         List<ALogPAtom> matchList = atomTypesByAtomSym.get(atSym);
         if( matchList == null )
         {  System.err.printf("Missing atom type for %s in %s\n",
                              Atom.getAtomName(at), oechem.OECreateSmiString(mol));
            continue;
         }

         for( ALogPAtom aType : matchList)
         {  if( aType.matchesAtom(at, hAreExplicit) )
            {  aLogP += aType.getHydrophobicity();
               //System.err.printf("At Contrib %s %f\n",aType.getDescription(), aType.getHydrophobicity());
               atomTypeCounts[aType.getType()]++;
               if( validatedAssignment) atomWasAssigned[at.GetIntData(ALOGPAtIdxTag)] = true;
               break;
            }
         }
      }
      atIt.delete();
      return aLogP;
   }

   /** get contribution of hydrogens on this atom
    */
   private double getAtomHContrib(OEAtomBase at, boolean hAreExplicit)
   {  int totlHCount = at.GetTotalHCount();
      if( totlHCount == 0 ) return 0D;

      List<ALogPAtom> matchList = atomTypesByAtomSym.get("H");
      for( ALogPAtom aType : matchList)
      {  if( aType.matchesAtom(at, hAreExplicit) )
         {  atomTypeCounts[aType.getType()] += totlHCount;
            if( validatedAssignment ) hWasAssigned[at.GetIntData(ALOGPAtIdxTag)] = true;

            //System.err.printf("atomHContrib %s %f\n", aType.getDescription(), aType.getHydrophobicity());
            return totlHCount * aType.getHydrophobicity();
         }
      }

      return 0D;
   }


   /**
    * return the atom counts for the last molecule passt to {@link #computeALogP}.
    *
    * To get the name of the count at each position call {@link #getAtomTypeName}.
    */
   public int[] getAtomCounts()
   {  return atomTypeCounts.clone();
   }

   private void validateAssignments(OEMolBase mol)
   {  OEAtomBaseIter atIt = mol.GetAtoms();
      while( atIt.hasNext() )
      {  OEAtomBase at = atIt.next();
         if( at.GetAtomicNum() != 1 )
         {  if( ! atomWasAssigned[at.GetIdx()] )
            {  System.err.printf("Atom %s in %s was not assigned\n",
                        Atom.getAtomName(at), oechem.OECreateSmiString(mol));
            }

            if( at.GetTotalHCount() > 0 && ! hWasAssigned[at.GetIdx()] )
            {  System.err.printf("H on %s in %s was not assigned\n",
                     Atom.getAtomName(at), oechem.OECreateSmiString(mol));
            }
         }
      }
      atIt.delete();
   }


   /** assign idx tag to each atom so that we can map back from molecule copies
    */
   private static void assignIndexTag(OEMolBase mol)
   {  OEAtomBaseIter atIt = mol.GetAtoms();
      while( atIt.hasNext() )
      {  OEAtomBase at = atIt.next();
         at.SetIntData(ALOGPAtIdxTag, at.GetIdx());
      }
      atIt.delete();
   }


   /**
    * Parse the tab-delimited file with the SMARTS definitions
    *
    * @return  the highest number of the atomType found.
    *
    */
   private int readSMARTS( String file )
   {  int smartsCounter = 0;
      int maxOfficialType = 0;

      try
      {  BufferedReader reader;

         if( file == null || file.length() == 0 )
         {  InputStream strm = this.getClass().getResourceAsStream("aLogPFragments.txt");
            reader = new BufferedReader(new InputStreamReader(strm));
         } else
         {   reader = new BufferedReader( new FileReader( file ) );
         }

         String line;

         while( null != ( line = reader.readLine() ) )
         {  line = line.trim();
            if( line.length() == 0 ||
                line.startsWith( "#" ) || line.startsWith( "\"#" ) )
               continue;

            String[] fields = line.split( "\t" ); //#AlogPName Type   Central Atom   Match Order   H Explicit  SMARTS   Description   hydrophobicity

            if( fields.length < 7 )
               throw new Error("Line in fragment file has less than 6 columns\n" + line);

            int aLogPType   = Integer.parseInt(fields[1].trim().replaceAll( "\"", "" ));
            String atom     = fields[2].trim().replaceAll( "\"", "" );
            float matchOrder= Float.parseFloat(fields[3].trim().replaceAll( "\"", "" ));
            boolean needExplicitH = "y".equalsIgnoreCase(fields[4].trim().replaceAll( "\"", "" ));
            String smarts   = fields[5].trim().replaceAll( "\"", "" );
            String desc     = fields[6].trim().replaceAll( "\"", "" );
            double hydroPhob= Double.parseDouble(fields[7].trim().replaceAll( "\"", "" ));

            ALogPAtom aLogPAtomType
               = new ALogPAtom( aLogPType, atom, matchOrder, smarts, needExplicitH, desc, hydroPhob );

            allAtomTypes.add(aLogPAtomType);

            if( aLogPType > maxOfficialType ) maxOfficialType = aLogPType;

            smartsCounter++;
         }
      } catch (IOException e)
      {  throw new Error(e);
      }

      Collections.sort(allAtomTypes);     //sort by matchorder

      for( ALogPAtom a : allAtomTypes)
      {  List<ALogPAtom> atTypeByAtList = atomTypesByAtomSym.get(a.getAtomSymbol());
         if( atTypeByAtList == null )
         {  atTypeByAtList = new ArrayList<ALogPAtom>();
            atomTypesByAtomSym.put(a.getAtomSymbol(), atTypeByAtList);
         }
         atTypeByAtList.add(a);
      }

      System.err.printf( "Read %d atom types\n", smartsCounter );

      return maxOfficialType;
   }

}

class ALogPAtom implements Comparable<ALogPAtom>
{  private final int aLogPType;
   private final String centerAtomSymbol;
   private final OESubSearch subSearch;
   private final boolean needExplicitH;
   private final String description;
   private final double hydroPhob;
   @SuppressWarnings("unused")
   private final String smarts;
   private final float matchOrder;

   public ALogPAtom(int aLogPType, String atom, float matchOrder, String smarts,
            boolean needExplicitH, String desc, double hydroPhob)
   {  this.aLogPType = aLogPType;
      this.matchOrder = matchOrder;
      this.centerAtomSymbol = atom;
      this.smarts = smarts;
      this.subSearch = new OESubSearch(smarts);
      this.needExplicitH = needExplicitH;
      this.description = desc;
      this.hydroPhob = hydroPhob;

      if( ! subSearch.IsValid() )
         throw new Error(String.format("Invalid smarts (%d): %s\n", aLogPType, smarts));
   }

   public String getAtomSymbol()
   {  return centerAtomSymbol;
   }

//   public int getIndex()
//   {  return idx;
//   }
//
   public double getHydrophobicity()
   {  return hydroPhob;
   }

   public String getDescription()
   {  return description; }

   public boolean matchesAtom(OEAtomBase at, boolean hAreExplicit )
   {  if( hAreExplicit != needExplicitH ) return false;
//System.err.println(subSearch.AtomMatch(at) + " " + at.GetTotalHCount()+" " +at.GetImplicitHCount()+ " " + at.GetExplicitHCount()+" " + at.GetDegree());
      return subSearch.AtomMatch(at);
   }

   /** number of this atom type as in paper */
   public int getType()
   {  return aLogPType;
   }

   public void close()
   {  subSearch.delete();
   }

   @Override
   public int compareTo(ALogPAtom o)
   {  return Float.compare(this.matchOrder, o.matchOrder);
   }
}
