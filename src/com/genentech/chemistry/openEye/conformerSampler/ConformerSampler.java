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
package com.genentech.chemistry.openEye.conformerSampler;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

import openeye.oechem.*;

import com.aestel.math.CircularRandomNumberGenerator;
import com.aestel.math.MixedBaseNumberConverter;

public class ConformerSampler
{  /** to make bonds fixed use bd.setBoolData(ISFixedBondTag, true) */
   public static final int  ISFixedBondTag = oechem.OEGetTag("TorScanner.bIsFixed");

   private static final String DefaultTorPatternFilename = "HRotorPattern.txt";

   private final TorsionPattern[] patterns;

   public ConformerSampler(String patternFile)
   {  patterns = parsePatternFile(patternFile);
   }

   /**
    * Generate maxConfs conformers by evaluating all combinatorial possibliities
    * of torsion around torsion described in the pattern file.
    *
    * Bonds cna be excluded from the sampling by using bd.setBoolData(ISFixedBondTag, true).
    *
    * @param mol input 3d mol
    * @param maxConfs if you wnat to be exhaustive use big number, will never be
    *        than total number of possiblilies
    * @return OEMCMolBase with conformers in random order.
    */
   public OEMCMolBase createConformations(OEMolBase mol, long maxConfs)
   {  OEMCMolBase confMol = oechem.OENewMCMolBase(mol);
      oechem.OEAddExplicitHydrogens(confMol, false, true);
      OEConfBaseIter confIt = confMol.GetConfs();
      OEConfBase baseConf = confIt.next();
      confIt.delete();

      List<Rotor> rotors = getRotors(confMol);

      int[] scanCounts = new int[rotors.size()];
      for(int i=0; i<rotors.size(); i++)
         scanCounts[i] = rotors.get(i).matchedPattern.rotations.length;

      MixedBaseNumberConverter mbnc = new MixedBaseNumberConverter(scanCounts);
      // subtract one because input conformation is always returned first
      long maxMolCount = Math.min(maxConfs-1, mbnc.getMaxValue());
      CircularRandomNumberGenerator crand = new CircularRandomNumberGenerator(mbnc.getMaxValue());

      if (maxMolCount > 99)
      {  String msg = String.format("WARNING: Expanding around fixed torsion with max of %d conformers.", maxMolCount);
         System.err.println(msg);
      }

      for(int molCount=0; molCount<maxMolCount; molCount++)
      {  OEConfBase newConf = confMol.NewConf(baseConf);

         long molNum;
         while((molNum = crand.nextValue()) == 0) // zeroed mol is returned on 0
            ;

         // get array of indexes into the specific torsion increment for each
         // found rotor
         int[] torsionIdxForRotor = mbnc.getMixedBasedDigits(molNum);
         for(int tor=0; tor<torsionIdxForRotor.length; tor++)
         {  Rotor rotor = rotors.get(tor);
            OEAtomBase a1 = rotor.atoms[0];
            OEAtomBase a2 = rotor.atoms[1];
            OEAtomBase a3 = rotor.atoms[2];
            OEAtomBase a4 = rotor.atoms[3];
            double angleIncrement = rotor.matchedPattern.rotations[torsionIdxForRotor[tor]];
            double angle = rotor.initialAngle + angleIncrement;
            oechem.OESetTorsion(newConf, a1, a2, a3, a4, angle);
         }
      }

      return confMol;
   }

   private List<Rotor> getRotors(OEMCMolBase confMol)
   {  boolean[] bondFlipped = new boolean[confMol.GetMaxBondIdx()+1];
      List<Rotor> rotors = new ArrayList<Rotor>();

      for(TorsionPattern tPat : patterns)
      {  OEMatchBaseIter matIt = tPat.subSearch.Match(confMol);
         while(matIt.hasNext())
         {  OEMatchBase m = matIt.next();

            OEMatchPairBondIter bIt = m.GetBonds();
            bIt.next();
            assert bIt.hasNext();
            OEBondBase b = bIt.next().getTarget();
            bIt.delete();

            if( b.GetBoolData(ISFixedBondTag) )
            {
               continue; // marged as fixed
            }
            if( bondFlipped[b.GetIdx()] ) continue; // bond already matched
            bondFlipped[b.GetIdx()] = true;

            OEMatchPairAtomIter atIt = m.GetAtoms();
            OEAtomBase[] torAtoms = new OEAtomBase[4];
            for(int i=0; i<4; i++)
            {  OEAtomBase at = atIt.next().getTarget();
               torAtoms[3-i] = at; // reverse so that first is last
            // the reason we want to have reverse order is that OESetTorsion
            // moves the last atom and we want to define the smarts such that
            // the first atom is moved
            }
            atIt.delete();

            double angle = oechem.OEGetTorsion(confMol, torAtoms[0], torAtoms[1],
                                                        torAtoms[2], torAtoms[3]);
            rotors.add(new Rotor(tPat, torAtoms, angle));
         }
         matIt.delete();
      }
      return rotors;
   }

   public void close()
   {  for( TorsionPattern p : patterns)
         p.close();
   }


   private TorsionPattern[] parsePatternFile(String patternFile)
   {  String line = null;
      int lineNum = 0;

      ArrayList<TorsionPattern> patList = new ArrayList<TorsionPattern>();
      try
      {  BufferedReader reader = null;

         if( patternFile != null && patternFile.length() > 0 )
         {
            // If the over-riding file exists, load it
            File f = new File(patternFile);
            if(f.exists() && !f.isDirectory())
            {
               System.err.println("ConformerSampler: Loading pattern file: " + patternFile);
               reader = new BufferedReader( new FileReader( patternFile ) );
            }
            else
            {
               // Else try to load the over-riding file as a resource in the class package
               InputStream strm = this.getClass().getResourceAsStream(patternFile);
               if (strm != null)
               {
                  System.err.println("ConformerSampler: Loading pattern file: " + patternFile);
                  reader = new BufferedReader(new InputStreamReader(strm));
               }
               else
               {  System.err.println("ConformerSampler: Failed to load pattern file: " + patternFile); }
            }
         }

         // Load default pattern if either patternFile == null and or the patternFile failed to load
         if (reader == null)
         {
            System.err.println("ConformerSampler: Loading default pattern file.");
            InputStream strm = this.getClass().getResourceAsStream(DefaultTorPatternFilename);
            reader = new BufferedReader(new InputStreamReader(strm));
         }

         while( null != ( line = reader.readLine() ) )
         {  lineNum++;
            line = line.trim();
            if( line.length() == 0 ||
                line.startsWith( "#" ) || line.startsWith( "\"#" ) )
               continue;

            String[] fields = line.split( "\t" ); //#Smarts rotations

            if( fields.length < 2 )
               throw new Error("Line in pattern file has less than 2 columns\n" + line);

            String smarts   = fields[0].trim().replaceAll( "\"", "" );
            String[] rotStr= fields[1].trim().replaceAll( "\"", "" ).split(",");
            double[] rot = new double[rotStr.length];
            for( int i=0; i<rotStr.length; i++)
               rot[i] = Math.toRadians(Double.parseDouble(rotStr[i]));

            patList.add(new TorsionPattern(smarts, rot));
         }

         reader.close();

         return patList.toArray(new TorsionPattern[patList.size()]);

      } catch(Throwable e)
      {  String msg = String.format("FileName=%s, line=%d error=%s",
               patternFile, lineNum, e.getMessage());
         throw new Error(msg, e);
      }
   }
}

class TorsionPattern
{  final OESubSearch subSearch;
   final String smarts;
   final double[] rotations;

   public TorsionPattern(String smarts, double[] rot)
   {  this.smarts = smarts;
      this.subSearch = new OESubSearch(smarts);
      this.rotations = rot.clone();
   }

   public void close()
   {  subSearch.delete();
   }
}

class Rotor
{  final TorsionPattern matchedPattern;
   final OEAtomBase[] atoms;
   final double initialAngle;

   Rotor(TorsionPattern tors, OEAtomBase[] atoms, double angle)
   {  assert atoms.length == 4;

      this.matchedPattern = tors;
      this.atoms = atoms.clone();
      this.initialAngle = angle;
   }
}
