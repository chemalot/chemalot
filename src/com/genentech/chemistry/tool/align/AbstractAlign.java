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
package com.genentech.chemistry.tool.align;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;

import com.aestel.utility.DataFormat;
import com.genentech.chemistry.tool.align.SDFAlign.OUTType;

import openeye.oechem.OEAtomBase;
import openeye.oechem.OEAtomBaseIter;
import openeye.oechem.OEElemNo;
import openeye.oechem.OEGraphMol;
import openeye.oechem.OEIsChiralAtom;
import openeye.oechem.OEIsHeavy;
import openeye.oechem.OEMatch;
import openeye.oechem.OEMatchBase;
import openeye.oechem.OEMatchBaseIter;
import openeye.oechem.OEMatchPairAtom;
import openeye.oechem.OEMatchPairAtomIter;
import openeye.oechem.OEMol;
import openeye.oechem.OEMolBase;
import openeye.oechem.OEUnaryAtomPred;
import openeye.oechem.oechem;

class AlignmentResult
{  public final double rmsd;
   public final int nMatch;
   public final int nHeavyRef;
   public final double[] atomDists;
   public final HashMap<String, String[]> atomInfo;
   
   public AlignmentResult(double rmsd, int nMatch, int nHeavyRef, double[] atomDists, HashMap<String, String[]> atomInfo)
   {  this.rmsd = rmsd;
      this.nMatch = nMatch;
      this.nHeavyRef = nHeavyRef;
      this.atomDists = atomDists;
      this.atomInfo  = atomInfo;
   }
}

abstract class AbstractAlign implements AlignInterface
{  public static final String ALIGN_SIZE_TAG = "AlignSize";
   private static final String ALIGN_PCT_TAG = "AlignPct";
   private static final OEIsChiralAtom CAF = new OEIsChiralAtom();
   protected final boolean doMirror;
   protected final String rmsdTag;
   private final String atomDevTag;
   private final OUTType oMolType;
   private final OEMolBase mCopy = new OEMol();
   protected final boolean doOptimize;
   private final boolean quiet;
   private String[] atomInfoTags;
   private final OEUnaryAtomPred IS_HeavyPred = new OEIsHeavy();

   @Override
   public abstract void close();

   public AbstractAlign(OUTType oMolType, String rmsdTag, String atomDevTag, String[] atomInfoTags, boolean doOptimize, boolean doMirror, boolean quiet)
   {  this.doMirror = doMirror;
      this.rmsdTag = rmsdTag;
      this.atomDevTag = atomDevTag;
      this.atomInfoTags = atomInfoTags;
      this.oMolType = oMolType;
      this.doOptimize = doOptimize;
      this.quiet = quiet;
   }


   public static void createMirror(OEMolBase fitmol)
   {  float coords[] = new float[fitmol.GetMaxAtomIdx() * 3];
      fitmol.GetCoords(coords);
      for( int i=0; i < coords.length; i+=3 )
         coords[i] = -coords[i];
      fitmol.SetCoords(coords);
   }


   public static boolean isChiral( OEMolBase mol )
   {  OEAtomBaseIter it = mol.GetAtoms(CAF);
      while(it.hasNext())
      {  OEAtomBase at = it.next();
         if( at.GetAtomicNum() == OEElemNo.N && at.GetDegree() != 4 )
            continue; // ignore stereo on N

         it.delete();
         return true;
      }
      it.delete();
      return false;
   }

   @Override
   public void align(OEMolBase fitmol)
   {  mCopy.Clear();
      oechem.OEAddMols( mCopy, fitmol );
      mCopy.SetTitle(fitmol.GetTitle());

      AlignmentResult bestAlign = doAlign(fitmol, Double.MAX_VALUE, true);

      oechem.OEPerceiveChiral(fitmol);
      if( doMirror && ! isChiral(fitmol) )
      {  OEGraphMol mirMol = new OEGraphMol(fitmol);
         createMirror(mirMol);

         AlignmentResult newAlign = doAlign(mirMol, bestAlign.rmsd, false);

         if( newAlign.rmsd < bestAlign.rmsd )
         {  fitmol.Clear();
            oechem.OEAddMols(fitmol, mirMol);
            fitmol.SetTitle(mirMol.GetTitle());
            bestAlign = newAlign;
         }
      }

      if( oMolType == OUTType.ORIGINAL )
      {  fitmol.Clear();
         oechem.OEAddMols(fitmol, mCopy);
         fitmol.SetTitle(mCopy.GetTitle());
      }

      if( bestAlign.rmsd < Double.MAX_VALUE && rmsdTag != null )
      {  oechem.OESetSDData(fitmol, rmsdTag, DataFormat.formatNumber(bestAlign.rmsd, "si3"));
         oechem.OESetSDData(fitmol, AbstractAlign.ALIGN_SIZE_TAG, Integer.toString(bestAlign.nMatch));
         
         //bestAlign.
         
         oechem.OESetSDData(fitmol, AbstractAlign.ALIGN_PCT_TAG, 
               DataFormat.formatNumber((100.*bestAlign.nMatch)/bestAlign.nHeavyRef, "r1"));
      }
      
      addAtomDistances(fitmol, bestAlign);
      addAtomMappedData(fitmol, bestAlign);
   }

   
   protected void addAtomDistances(OEMolBase fitmol, AlignmentResult bestAlign)
   {  if( bestAlign.atomDists != null )
      {  StringBuilder sb = new StringBuilder();
         for( double d:  bestAlign.atomDists)
         {  String dStr = d < 0 ? "" : DataFormat.formatNumber(d, "si3");  
            sb.append(dStr).append(',');
         }
         if( sb.length() > 0 ) sb.setLength(sb.length()-1);  //last ","
         
         String distStr = sb.toString()
                            .replaceAll("\\.?0+,", ",");
                            
         oechem.OESetSDData(fitmol, atomDevTag, distStr);
      }
   }

   
   protected void addAtomMappedData(OEMolBase fitmol, AlignmentResult bestAlign)
   {  if( bestAlign.atomInfo != null)
      {  for( Entry<String, String[]> aie : bestAlign.atomInfo.entrySet() )
         {  String tag = aie.getKey();
            String[] vals = aie.getValue();
            
            boolean foundVal = false;
            StringBuilder sb = new StringBuilder();
            for( String v : vals)
            {  if( v != null && v.length() > 0) foundVal = true;
               sb.append(v == null ? "" : v).append(",");
            }
            
            if(foundVal)
            {  if(sb.length() > 1) sb.setLength(sb.length()-1);
               oechem.OESetSDData(fitmol, tag, sb.toString());
            }
         }
      }
   }

   
   /** list of reference molecules to which to align.
    *  This method must return the same number of elements as {@link #getReference()}
    *  and the elements must be in the same order
    **/
   abstract protected List<OEMolBase> getReference();

   /** List of MatchIterators that match Reference to fitMol.
    *  This method must return the same number of elements as {@link #getReference()}
    *  and the elements must be in the same order
    **/
   abstract public List<OEMatchBaseIter> getMatchIterator(OEMolBase fitMol);
   

   private AlignmentResult doAlign(OEMolBase fitMol, double minRMSD, boolean warnOnNoMatch)
   {  double[] bestRMat = null;
      double[] bestTrans = null;
      OEMolBase bestRef = null;
      OEMatchBase bestMatch = null;
      int nHeavyRef = -1;
      int nMatchAtMin = -1;
      
      
      Iterator<OEMolBase> refs = getReference().iterator();
      for( OEMatchBaseIter matchIt : getMatchIterator(fitMol))
      {  OEMolBase ref = refs.next();

         while( matchIt.hasNext() )
         {  OEMatchBase match = matchIt.next();

            // we might have hydrogens in the core but want to calculate the rmsd
            // only on the heavy atoms
            OEMatch nonHMatch = getNonHMatch(match);
            match.delete();
            if( nonHMatch.NumAtoms() == 0 )
            {  nonHMatch.delete();
               continue;
            }

            double rmat[] = null;
            double trans[] = null;
            double rmsd;
            if( doOptimize )
            {  rmat = new double[9];
               trans = new double[3];
               rmsd = oechem.OERMSD(ref, fitMol, nonHMatch, doOptimize, rmat, trans);
            } else
            {  rmsd = oechem.OERMSD(ref, fitMol, nonHMatch, doOptimize);
            }
            
            if( rmsd < minRMSD )
            {  bestRMat = rmat;
               bestTrans = trans;
               minRMSD = rmsd;
               nMatchAtMin = nonHMatch.NumAtoms();
               nHeavyRef = oechem.OECount(ref, IS_HeavyPred);
               
               if( atomDevTag != null || atomInfoTags.length > 0)
               {  if( bestMatch != null ) 
                  {  bestMatch.delete();
                     if( bestRef != ref ) bestRef.delete();
                  }
                  bestMatch = nonHMatch;
                  bestRef = ref;
               } else
               {  nonHMatch.delete();
               }
            } else
            {  nonHMatch.delete();
            }

         }
         if( bestRef != ref ) ref.delete();
         matchIt.delete();
      }

      if( bestRMat != null && doOptimize )
      {  oechem.OERotate(fitMol, bestRMat);
         oechem.OETranslate(fitMol, bestTrans);
      } else
      {  if(! quiet && nMatchAtMin == -1 && warnOnNoMatch) System.err.println("sdfAlign: No match");
      }
      
      double[] atomDists = null;
      HashMap<String,String[]>atomInfo = null;
      
      if( bestMatch != null )
      {  atomDists = atomDevTag == null ? null : new double[fitMol.NumAtoms()];
         if( atomDevTag != null) Arrays.fill(atomDists, -1.);
         
         atomInfo = atomInfoTags.length > 0 ? new HashMap<String,String[]>() : null;
            
         OEMatchPairAtomIter matchAtIt = bestMatch.GetAtoms();
         while( matchAtIt.hasNext() )
         {  OEMatchPairAtom matchPair = matchAtIt.next();
            OEAtomBase refAt = matchPair.getPattern();
            OEAtomBase molAt = matchPair.getTarget();
            
            if( atomDevTag != null) 
               atomDists[molAt.GetIdx()] = oechem.OEGetDistance(bestRef, refAt, fitMol, molAt);
            
            if( atomInfoTags.length > 0 )
            {  for( String tag : atomInfoTags )
               {  String atData = refAt.GetStringData(tag);
                  if( atData != null )
                  {  String[] mapData = atomInfo.get(tag);
                     if(mapData == null)
                     {  mapData = new String[fitMol.NumAtoms()];
                        atomInfo.put(tag, mapData);
                     }
                     mapData[molAt.GetIdx()] = atData;
                  }
               }
            }
         }
         bestMatch.delete();
         bestRef.delete();
      }

      return new AlignmentResult(minRMSD, nMatchAtMin, nHeavyRef, atomDists, atomInfo);
   }

   private static OEMatch getNonHMatch(OEMatchBase match)
   {  OEMatch nonHMatch = new OEMatch();
      OEMatchPairAtomIter atIt = match.GetAtoms();

      while(atIt.hasNext())
      {  OEMatchPairAtom at = atIt.next();
         if(at.getTarget().GetAtomicNum() != 1)
         {  nonHMatch.AddPair(at.getPattern(), at.getTarget());
         }
      }
      atIt.delete();
      return nonHMatch;
   }
}

