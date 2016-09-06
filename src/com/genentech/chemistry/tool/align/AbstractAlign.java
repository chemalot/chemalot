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

import java.util.Iterator;
import java.util.List;

import openeye.oechem.*;

import com.aestel.utility.DataFormat;
import com.genentech.chemistry.tool.align.SDFAlign.OUTType;

abstract class AbstractAlign implements AlignInterface
{  static OEIsChiralAtom CAF = new OEIsChiralAtom();
   protected final boolean doMirror;
   protected final String rmsdTag;
   private final OUTType oMolType;
   private final OEMolBase mCopy = new OEMol();
   protected final boolean doOptimize;
   private final boolean quiet;

   @Override
   public abstract void close();

   public AbstractAlign(OUTType oMolType, String rmsdTag, boolean doOptimize, boolean doMirror, boolean quiet)
   {  this.doMirror = doMirror;
      this.rmsdTag = rmsdTag;
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

      double minRMSD = doAlign(fitmol, Double.MAX_VALUE, true);

      oechem.OEPerceiveChiral(fitmol);
      if( doMirror && ! isChiral(fitmol) )
      {  OEGraphMol mirMol = new OEGraphMol(fitmol);
         createMirror(mirMol);

         double newRMSD = doAlign(mirMol, minRMSD, false);

         if( newRMSD < minRMSD )
         {  fitmol.Clear();
            oechem.OEAddMols(fitmol, mirMol);
            fitmol.SetTitle(mirMol.GetTitle());
            minRMSD = newRMSD;
         }
      }

      if( oMolType == OUTType.ORIGINAL )
      {  fitmol.Clear();
         oechem.OEAddMols(fitmol, mCopy);
         fitmol.SetTitle(mCopy.GetTitle());
      }

      if( minRMSD < Double.MAX_VALUE && rmsdTag != null )
            oechem.OESetSDData(fitmol, rmsdTag, DataFormat.formatNumber(minRMSD, "si3"));

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

   private double doAlign(OEMolBase fitMol, double minRMSD, boolean warnOnNoMatch)
   {  double[] bestRMat = null;
      double[] bestTrans = null;

      Iterator<OEMolBase> refs = getReference().iterator();
      for( OEMatchBaseIter matchIt : getMatchIterator(fitMol))
      {  OEMolBase ref = refs.next();

         while( matchIt.hasNext() )
         {  OEMatchBase match = matchIt.next();

            // we might have hydrogens in the core but want to calculate the rmsd
            // only on the heavy atoms
            OEMatch nonHMatch = getNonHMatch(match);

            double rmat[] = new double[9];
            double trans[] = new double[3];


            double rmsd = oechem.OERMSD(ref, fitMol, nonHMatch, doOptimize, rmat, trans);
            if( rmsd < minRMSD )
            {  bestRMat = rmat;
               bestTrans = trans;
               minRMSD = rmsd;
            }

            nonHMatch.delete();
            match.delete();
         }
         ref.delete();
         matchIt.delete();
      }

      if( bestRMat != null && doOptimize )
      {  oechem.OERotate(fitMol, bestRMat);
         oechem.OETranslate(fitMol, bestTrans);
      } else
      {  if(! quiet && warnOnNoMatch) System.err.println("sdfAlign: No match");
      }
      return minRMSD;
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

