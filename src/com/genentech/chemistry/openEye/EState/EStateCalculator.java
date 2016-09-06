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
package com.genentech.chemistry.openEye.EState;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import openeye.oechem.OEAtomBase;
import openeye.oechem.OEAtomBaseIter;
import openeye.oechem.OEBondBase;
import openeye.oechem.OEBondBaseIter;
import openeye.oechem.OEMatchBase;
import openeye.oechem.OEMatchBaseIter;
import openeye.oechem.OEMolBase;
import openeye.oechem.OESubSearch;
import openeye.oechem.oechem;

import com.aestel.chemistry.molecule.Atom;
import com.aestel.utility.NameValuePair;



/**
 * Compute the EState value and atom group assignment based on
 * 
 * Hall, L. H.; Kier, L. B. Electrotopological State Indices for Atom Types:
 * A Novel Combination of Electronic, Topological, and Valence State Information.
 * J. Chem. Inf. Comput. Sci. 1995, 35, 1039-1045
 *
 * @author Man-Ling Lee / August 2, 2012
 * Copyright 2012-2015 Genentech
 */
public class EStateCalculator
{
   OEMolBase mol;
   List<EStateCalcResult> esCalcResultList;
   EStateCalcResult[] esCalcResults;
   float[] estateAtomGroupSums;
   int[] estateAtomGroupCounts;
   int notAssignedAtomCounter;
   String notAssignedAtomSummary;
   
   
   public EStateCalculator()
   {  esCalcResultList      = new ArrayList<EStateCalcResult>();
      estateAtomGroupSums   = new float[ EStateAtomGroup.getEStateAtomGroupSize() ];
      estateAtomGroupCounts = new int[ EStateAtomGroup.getEStateAtomGroupSize() ];
      clear();
   }
   
   
   private void clear()
   {  esCalcResultList.clear();
      mol = null;
      for( int i=0; i<estateAtomGroupSums.length; i++ )
         estateAtomGroupSums[i] = 0.0f;
      for( int i=0; i<estateAtomGroupCounts.length; i++ )
         estateAtomGroupCounts[i] = 0;
      notAssignedAtomCounter = 0;
      notAssignedAtomSummary = "";
   }
   
   
   public void compute( OEMolBase mol, boolean doKeys, boolean printDetails )
   {  this.clear();
      this.mol = mol;
      assignEStateAtomGroup( printDetails );
      //if( doKeys )
         computeEStateSums( printDetails );
   }
   
   
   public String getEStateIndexSummary()
   {  StringBuilder sb = new StringBuilder();
      for( EStateCalcResult result : esCalcResultList )
      {  sb.append( result.getEStateIndex() ).append( " (" )
           .append( result.getAtomLabel() ).append( "); " );
      }
      return sb.substring( 0, sb.length()-2 );
   }
   
   
   @SuppressWarnings("unchecked")
   public NameValuePair<String,String>[] getEStateOf( OEMolBase mol, String smarts )
   {  OESubSearch subSearch = new OESubSearch( smarts );
      OEMatchBaseIter matchIt = subSearch.Match( mol, true );
      
      List<NameValuePair<String,String>> nvpList = new ArrayList<NameValuePair<String,String>>();
      if( matchIt.hasNext() )
      {  NameValuePair<String,String> nvp;
         OEMatchBase match = matchIt.next();
         OEAtomBaseIter targetIt  = match.GetTargetAtoms();
         OEAtomBaseIter patternIt = match.GetPatternAtoms();
         while( targetIt.hasNext() && patternIt.hasNext() )
         {
            OEAtomBase tAtom = targetIt.next();
            OEAtomBase pAtom = patternIt.next();
            EStateCalcResult result = esCalcResults[ tAtom.GetIdx() ];
            String estate = String.valueOf( result.getEStateIndex() );
            
            String tag = "Index_" + ( pAtom.GetIdx()+1 );
            nvp = new NameValuePair<String,String>( tag, estate );
            nvpList.add( nvp );
            
            tag = tag + "_atom";
            nvp = new NameValuePair<String,String>( tag, result.getAtomLabel() );
            nvpList.add( nvp );
            
            pAtom.delete();
            tAtom.delete();
         }
         targetIt.delete();
         patternIt.delete();
         match.delete();
      }
      matchIt.delete();
      subSearch.delete();
      return nvpList.toArray( new NameValuePair[ nvpList.size() ] );
   }
   
   
   /**
    * Returns an array containing the occurrences for each E-state atom group.
    * Note that it a copy.
    */
   public int[] getEStateCounts()
   {  return estateAtomGroupCounts.clone(); }
   
   /**
    * Returns an array containing the sums for each E-state atom group.
    * Note that it a copy.
    */
   public float[] getEStateSums()
   {  return estateAtomGroupSums.clone(); 
   }
   
   /**
    * Returns an array containing the E-state atom group symbols.
    * Note that it a copy.
    */
   public String[] getEStateAtomGroupSymbols()
   {  return EStateAtomGroup.getEStateAtomGroupSymbols(); }
   
   /**
    * Returns the number of atoms in the given molecules could not be assigned
    * to an atom group.
    */
   public int getUnknownCount()
   {  return notAssignedAtomCounter; }


   public String getUnknownAtoms()
   {  return this.notAssignedAtomSummary; }

   
   private void assignEStateAtomGroup( boolean printDetails )
   {  OEAtomBaseIter atomIterator = mol.GetAtoms();
      esCalcResults = new EStateCalcResult[ mol.GetMaxAtomIdx() ];
      
      while( atomIterator.hasNext() )
      {  EStateAtomGroup esAtomGroup = null;
         OEAtomBase atom = atomIterator.next();
         int atomicNum = atom.GetAtomicNum();
         if( atomicNum == 1 )
            continue;
         
         EStateCalcResult esCalcResult = new EStateCalcResult( atom, printDetails );
         esCalcResultList.add( esCalcResult );
         esCalcResults[ atom.GetIdx() ] = esCalcResult;
         
         int hash = EStateAtomGroup.computeHashNumber( 
                  esCalcResult.getAtomicNumber(), esCalcResult.getDeltaValence(), 
                  esCalcResult.getDeltaSimple(), esCalcResult.getSumOfDelta(), 
                  esCalcResult.getDifferenceOfDelta(), esCalcResult.isAromatic() );
         List<EStateAtomGroup> esAtomGroupList;
         esAtomGroupList= EStateAtomGroup.getEStateAtomGroups( hash );
         if( esAtomGroupList == null )
         {  ++notAssignedAtomCounter;
            String atomSymbol = oechem.OEGetAtomicSymbol( atomicNum );
            int seqNo = atom.GetIdx() + 1;
            String oldSummary = notAssignedAtomSummary;
            notAssignedAtomSummary = oldSummary + atomSymbol + seqNo + ";";
         
         } else
         {  if( esAtomGroupList.size() == 1 )
            {  esAtomGroup = esAtomGroupList.get( 0 );
               esCalcResult.setEStateAtomGroup( esAtomGroup );
               int esAtomGroupIndex = esAtomGroup.atomGroupNumber-1;
               ++estateAtomGroupCounts[ esAtomGroupIndex ];
            
            } else
            {  int esvSingleBondCount = esCalcResult.computeNumberOfSingleBonds();
               for( int i=0; i<esAtomGroupList.size(); i++ )
               {  esAtomGroup = esAtomGroupList.get( i );
                  if( esAtomGroup.countSingleBonds && 
                      esAtomGroup.singleBonds == esvSingleBondCount )
                  {  esCalcResult.setEStateAtomGroup( esAtomGroup );
                     int esAtomGroupIndex = esAtomGroup.atomGroupNumber-1;
                     ++estateAtomGroupCounts[ esAtomGroupIndex ];
                     break;
                  }
               }
            }
         }
         if( printDetails && esCalcResult.getEStateAtomGroup() != null )
         {  EStateAtomGroup esag = esCalcResult.getEStateAtomGroup();
            System.err.println( "ES Atom Group Num___" + esag.atomGroupNumber );
            System.err.println( "ES Atom Group Sym___" + esag.atomGroupSymbol );
         }
      }
      atomIterator.delete();
   }
   
   
   private void computeEStateSums( boolean printDetails )
   {  for( EStateCalcResult esCalcResult : esCalcResultList )
      {  String atomSymbol = oechem.OEGetAtomicSymbol( esCalcResult.getAtomicNumber() );
         int seqNo = esCalcResult.getAtomIndex() + 1;
         if( printDetails )
            System.err.println( "\nAtom " + atomSymbol + seqNo );
         
         esCalcResult.computeEStateIndex( esCalcResultList, printDetails );
         if( esCalcResult.getEStateAtomGroup() == null )
            continue;
         EStateAtomGroup esag = esCalcResult.getEStateAtomGroup();
         int esAGIndex = esag.atomGroupNumber-1;
         float oldESSum = estateAtomGroupSums[ esAGIndex ];
         float myESIndex = esCalcResult.getEStateIndex();
         estateAtomGroupSums[ esAGIndex ] = oldESSum + myESIndex;
         
         if( printDetails )
         {  String esagSymbol = esag.atomGroupSymbol;
            System.err.println( "Atom Group Sum (" + esagSymbol + ")___" + 
                                estateAtomGroupSums[ esAGIndex ] );
         }
      }
   }
}


class EStateAtomGroup
{  private static final EStateAtomGroup[] ATOM_GROUP_LIST;
   private static final Map<Integer,List<EStateAtomGroup>> ATOM_GROUP_MAP;
   static 
   {  String estateFile = "EStateAtomGroups.txt";
      Map<Integer,List<EStateAtomGroup>> esagMap = new HashMap<Integer,List<EStateAtomGroup>>();
      List<EStateAtomGroup> esagList = new ArrayList<EStateAtomGroup>();
      try
      {  InputStream strm = EStateAtomGroup.class.getResourceAsStream( estateFile );
         BufferedReader reader = new BufferedReader( new InputStreamReader( strm ) );
         String line;
         while( null != ( line = reader.readLine() ) )
         {  line = line.trim();
            if( line.length() == 0 || 
                line.startsWith( "#" ) || line.startsWith( "\"#" ) )
               continue;
            
            String[] fields = line.split( "\t" ); 
            //no., atom group, Z, dv, d, dv+d, dv-d, AR, count single bond, single bonds, group symbol, 
            //Hash number (for tests)
            if( fields.length < 11 )
               throw new Error("Line in estate file has less than 11 columns\n" + line);
   
            int groupNumber      = Integer.parseInt( fields[0].trim().replaceAll( "\"", "" ) );
            String groupName    = fields[1].trim().replaceAll( "\"", "" );
            int atomicNumber    = Integer.parseInt( fields[2].trim().replaceAll( "\"", "" ) );
            int deltaV          = Integer.parseInt(fields[3].trim().replaceAll( "\"", "" ) );
            int delta           = Integer.parseInt( fields[4].trim().replaceAll( "\"", "" ) );
            int sumOfDelta      = Integer.parseInt( fields[5].trim().replaceAll( "\"", "" ) );
            int diffOfDelta     = Integer.parseInt( fields[6].trim().replaceAll( "\"", "" ) );
            boolean isAromatic  = "1".equals( fields[7].trim().replaceAll( "\"", "") );
            boolean countSBonds = "y".equalsIgnoreCase( fields[8].trim().replaceAll( "\"", "" ) );
            int singleBonds = -1;
            String s = fields[9].trim().replaceAll( "\"", "" );
            if( s.length() > 0 )
               singleBonds = Integer.parseInt( s );
            String groupSymbol  = fields[10].trim().replaceAll( "\"", "" );
            
            EStateAtomGroup esag = new EStateAtomGroup( groupName, groupSymbol, 
                     groupNumber, atomicNumber, deltaV, delta, sumOfDelta, diffOfDelta,
                     isAromatic, countSBonds, singleBonds );
            esagList.add( esag );
            
            int hashNumber = esag.hashNumber;
            if( !esagMap.containsKey( hashNumber ) )
               esagMap.put( hashNumber, new ArrayList<EStateAtomGroup>() );
            esagMap.get( hashNumber ).add( esag );
         }
      } catch (IOException e)
      {  throw new Error(e);
      }
      ATOM_GROUP_MAP   = esagMap;
      ATOM_GROUP_LIST  = esagList.toArray( new EStateAtomGroup[ esagList.size() ] );
   }
   
   
   static int getEStateAtomGroupSize()
   {  return ATOM_GROUP_LIST.length; }
   
   
   static String[] getEStateAtomGroupSymbols()
   {  String[] symbols = new String[ ATOM_GROUP_LIST.length ];
      for( int i=0; i<symbols.length; i++ )
         symbols[i] = ATOM_GROUP_LIST[i].atomGroupSymbol;
      return symbols;
   }
   
   
   static List<EStateAtomGroup> getEStateAtomGroups( int hashNumber )
   {  return ATOM_GROUP_MAP.get( hashNumber ); }
   
   static int computeHashNumber( int atomicNumber, int deltaV, int delta,
            int sumOfDelta, int diffOfDelta, boolean isAromatic )
   {  //hashNumber = (C3*1000000)+(D3*100000)+(E3*10000)+(F3*100)+(G3+2)*10+H3
      return ( atomicNumber * 1000000 ) + 
             ( deltaV * 100000 ) + 
             ( delta * 10000 ) +
             ( sumOfDelta * 100 ) + 
             ( ( diffOfDelta + 2 ) * 10 ) +
             ( isAromatic ? 1 : 0 );
   }   
   

   final String atomGroupName;
   final String atomGroupSymbol;
   final int atomGroupNumber;
   final int atomicNumber;
   final int deltaV;
   final int delta;
   final int sumOfDelta;
   final int diffOfDelta;
   final boolean isAromatic;
   final boolean countSingleBonds;
   final int singleBonds;
   final int hashNumber;
   
   
   EStateAtomGroup( String groupName, String groupSymbol, int groupNumber, 
                    int atomicNumber, int deltaV, int delta, 
                    int sumOfDelta, int diffOfDelta, boolean isAromatic, 
                    boolean countSingleBonds, int singleBonds )
   {  this.atomGroupNumber  = groupNumber;
      this.atomGroupName    = groupName;
      this.atomGroupSymbol  = groupSymbol;
      this.atomicNumber     = atomicNumber;
      this.deltaV           = deltaV;
      this.delta            = delta;
      this.sumOfDelta       = sumOfDelta;
      this.diffOfDelta      = diffOfDelta;
      this.isAromatic       = isAromatic;
      this.countSingleBonds = countSingleBonds;
      this.singleBonds      = singleBonds;
      this.hashNumber       = this.computeHashNumber();
   }
   
   private int computeHashNumber()
   {  return EStateAtomGroup.computeHashNumber( atomicNumber, 
               deltaV, delta, sumOfDelta, diffOfDelta, isAromatic );
   }   
}



class EStateCalcResult
{  private OEAtomBase atom = null;
   /*
    * Simple Connectivity Value: 
    * Single bonds minus H-bonds corresponds to number of heavy neighbor atom */
   private int delta;
   /*
    * Valence Connectivity Value: 
    * Sum of sigma and pi bonds, lone pairs minus H-bonds */
   private int deltaV;
   /*
    * Intrincis Atomic State:
    * Based on the Kier-Hall electronegativity and derived from the ratio of 
    * that electronegativity to the number of skeletal sigma bonds for a given
    * atom */
   private float intrinsicState; 
   private float eStateIndex;
   private EStateAtomGroup atomGroup = null;

   
   EStateCalcResult( OEAtomBase atom, boolean printDetails )
   {  if( atom == null )
         throw new NullPointerException();
      int atomicNum = atom.GetAtomicNum();
      int pqNumber = Atom.PERIOD[ atomicNum ];
   
      this.atom = atom; 
      delta  = atom.GetHvyDegree();
      deltaV = Atom.SP_ELECTRONS[ atomicNum ] - atom.GetTotalHCount();
      intrinsicState = (float)( ( Math.pow( 2D/pqNumber, 2D ) * deltaV + 1D ) 
                       / (double)delta );
      
      if( printDetails )
      {  String atomSymbol = oechem.OEGetAtomicSymbol( atom.GetAtomicNum() );
         int seqNo = atom.GetIdx() + 1;
         System.err.println( "\nAtom " + atomSymbol + seqNo );
         System.err.println( "deltaV___" + deltaV );
         System.err.println( "delta___" + delta );
         System.err.println( "N___" + pqNumber );
         System.err.println( "iState___" + intrinsicState +
                             " = ( (2/" + pqNumber + ")^2 * " + deltaV + 
                             " + 1 ) / " + delta );
      }
   }

   
   /**
    * @return this index of the given atom in the molecule.
    */
   int getAtomIndex()
   {  return atom.GetIdx(); }
   
   int getAtomicNumber()
   {  return atom.GetAtomicNum(); }
   
   /**
    * @return the atom symbol plus the internal index in, e.g. following format C_1
    */
   String getAtomLabel()
   {  String symbol =  oechem.OEGetAtomicSymbol( atom.GetAtomicNum() );
      return symbol + "_" + ( atom.GetIdx()+1 );
   }
   
   int getDeltaSimple()
   {  return delta; }
   
   int getDeltaValence()
   {  return deltaV; }
   
   int getSumOfDelta()
   {  return deltaV + delta; }
   
   int getDifferenceOfDelta()
   {  return deltaV - delta; }
   
   boolean isAromatic()
   {  return atom.IsAromatic(); }
   
   /**
    * The intrinsic state of an atom, i.e. I = ((2/N)^2 * deltaV + 1) / delta
    */
   float getIntrinsicState()
   {  return this.intrinsicState; }
   
   void setEStateAtomGroup( EStateAtomGroup atomGroup )
   {  this.atomGroup = atomGroup; }
   
   EStateAtomGroup getEStateAtomGroup()
   {  return this.atomGroup; }
   
   float getEStateIndex()
   {  return this.eStateIndex; }
   
   
   int computeNumberOfSingleBonds()
   {  int singleBondCount = 0;
      OEBondBaseIter bondIterator = atom.GetBonds();
      while( bondIterator.hasNext() )
      {
         OEBondBase bond = bondIterator.next();
         if( bond.GetOrder() == 1 && !bond.IsAromatic() )
            ++singleBondCount;
      }
      bondIterator.delete();
      return singleBondCount;
   }
   
   
   /**
    * Calculate the E-state index of the given atom, i.e. 
    *          S[i] = I[i] + sumOf( I[i] - I[j]) / r[ij]^2 )
    * @param esCalcResults List of atoms in the same molecules with intrinsic state
    *          for calculation of the E-state index of the given atom.
    */
   void computeEStateIndex( List<EStateCalcResult> esCalcResults,
                            boolean printDetails )
   {  float sumOfPertubation = 0;
      float myState = this.getIntrinsicState();
      for( EStateCalcResult other : esCalcResults )
      {  int rij = oechem.OEGetPathLength(this.atom, other.atom ) + 1;
         float oldSum = sumOfPertubation;
         float otherState = other.getIntrinsicState();
         sumOfPertubation = (float) (oldSum + ( myState - otherState ) / Math.pow( rij, 2 ));
         
         if( printDetails )
         {  System.err.println( "sumOfPertubation___" + sumOfPertubation + 
                     " = " + oldSum + " + (" + myState + " - " + otherState + 
                     ") / " + rij + "^2" );
         }
      }
      eStateIndex = myState + sumOfPertubation;
      if( printDetails )
         System.err.println( "Estate index___" + eStateIndex );
   }
}


