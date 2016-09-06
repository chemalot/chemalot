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
package com.genentech.chemistry.openEye;

import openeye.oechem.OEAtomBase;

/**
 * This HeadAtomComputer computes similarities of atom types that should be
 * intuitive to chemists.
 *
 * THe implementation is purley empirical making eg. halogens more similar t eachother.
 *
 * @author albertgo
 *
 */
public class FuzzyHeadAtomComputer extends HeadAtomComputer
{  @SuppressWarnings("hiding")
   public static final HeadAtomComputer INSTANCE = new FuzzyHeadAtomComputer();

   private static final int CARBON   = 6;
   private static final int NITROGEN = 7;
   private static final int OXYGEN   = 8;
   private static final int FLOURINE = 9;
   private static final int SULPHOR  = 16;
   private static final int CHLORINE = 17;
   private static final int BROMINE  = 35;


   @Override
   protected double getHeadAtomSim(int atIdx1, IAAPathComputerInterface m1,
            int atIdx2, IAAPathComputerInterface m2)
   {  int aType1 = m1.getAtomType(atIdx1);
      int aType2 = m2.getAtomType(atIdx2);
      int atNum1 = m1.getAtomNum(aType1);
      int atNum2 = m2.getAtomNum(aType2);

      // ensure atomic number of at1 is smaller= than at2
      // ensure symmetry and that we do not need to implement both directions
      if( atNum1 > atNum2 )
      {  int d=atNum1;
         atNum1 = atNum2;
         atNum2 = d;

         d=aType1;
         aType1 = aType2;
         aType2 = d;

         d=atIdx1;
         atIdx1 = atIdx2;
         atIdx2 = d;

         IAAPathComputerInterface md = m1;
         m1 = m2;
         m2 = md;
      }
      OEAtomBase at1 = m1.getAtom(atIdx1);
      OEAtomBase at2 = m2.getAtom(atIdx2);

      switch( atNum1 )
      {
         case CARBON:   return computeCarbon(  at1, at2, aType1, aType2, atNum2);
         case NITROGEN: return computeNitrogen(at1, at2, aType1, aType2, atNum2);
         case OXYGEN:   return computeOxygen(  at1, at2, aType1, aType2, atNum2);
         case FLOURINE: return computeFlourine(at1, at2, atNum2);
         case SULPHOR:  return computeSulphor( at1, at2, atNum2);
         case CHLORINE: return computeChlorine(at1, at2, atNum2);
         case BROMINE:  return computeBromine( at1, at2, atNum2);

         default:
            if(aType1 == aType2) return 1D;
            return at1.GetHvyDegree() == at2.GetHvyDegree() ? 0.2D : 0D;
      }
   }


   private static double computeCarbon(OEAtomBase at1, OEAtomBase at2, int aType1,
            int aType2, int atNum2)
   {  if( atNum2 != CARBON )
      {  if( IAAPathGenerator.isAromatic(aType1) == IAAPathGenerator.isAromatic(aType2))
            return 0.3D;

         return at1.GetHvyDegree() == at2.GetHvyDegree() ? 0.2D : 0D;
      }

      /////// CARBON - CARBON
      if( aType1 == aType2 )
         return 1D;

      return 0.8D;
   }


   private static double computeNitrogen(OEAtomBase at1, OEAtomBase at2, int aType1,
            int aType2, int atNum2)
   {  if( atNum2 != NITROGEN )
      {  if( atNum2 == OXYGEN )
         {  if( IAAPathGenerator.isAromatic(aType1) == IAAPathGenerator.isAromatic(aType2) )
               return 0.5D;

            return 0.25D;
         }

         // non N non O
         if( IAAPathGenerator.isAromatic(aType1) == IAAPathGenerator.isAromatic(aType2))
            return 0.3D;

         return at1.GetHvyDegree() == at2.GetHvyDegree() ? 0.2D : 0D;
      }

      /////// Nitrogen Nitrogen
      if( aType1 == aType2 && at1.GetValence() == at2.GetValence() )
         return 1D;

      return 0.6D;
   }


   private static double computeOxygen(OEAtomBase at1, OEAtomBase at2, int aType1,
            int aType2, int atNum2)
   {  if( atNum2 != OXYGEN )
      {  if( atNum2 == SULPHOR )
         {  if( IAAPathGenerator.isAromatic(aType1) == IAAPathGenerator.isAromatic(aType2)
                && at1.GetValence() == at2.GetValence() )
               return 0.5D;

            return 0.0D;
         }

         if( IAAPathGenerator.isAromatic(aType1) == IAAPathGenerator.isAromatic(aType2))
            return  0.3D;

         return at1.GetHvyDegree() == at2.GetHvyDegree() ? 0.2D : 0D;
      }


      /////// OXYGEN - OXYGEN
      if( aType1 == aType2 )
         return 1D;

      return 0.8;
   }


   private static double computeFlourine(OEAtomBase at1, OEAtomBase at2, int atNum2)
   {
      if( atNum2 != FLOURINE )
      {  if( atNum2 == CHLORINE && at1.GetValence() == at2.GetValence() )
            return 0.5D;

         return at1.GetHvyDegree() == at2.GetHvyDegree() ? 0.2D : 0D;
      }

      /////// FLOURINE FLOURINE
      return 1D;
   }


   private static double computeSulphor(OEAtomBase at1, OEAtomBase at2, int atNum2)
   {
      if( atNum2 != SULPHOR )
         return at1.GetHvyDegree() == at2.GetHvyDegree() ? 0.2D : 0D;

      int v1 = at1.GetValence();
      int v2 = at2.GetValence();
      if( v1 == v2 ) return 1D;
      if( v1 > 2 && v2 > 2 ) return 0.8;

      return 0D;
   }


   private static double computeChlorine(OEAtomBase at1, OEAtomBase at2, int atNum2)
   {
      if( atNum2 != CHLORINE )
      {  if( atNum2 == BROMINE )
         {  if( at1.GetValence() == at2.GetValence() )
               return 0.6D;
            else
               return 0D;
         }

         return at1.GetHvyDegree() == at2.GetHvyDegree() ? 0.2D : 0D;
      }

      /////// CHLORINE CHLORINE
      if( at1.GetValence() == at2.GetValence() )
         return 1D;

      return 0.0D;
   }


   private static double computeBromine(OEAtomBase at1, OEAtomBase at2, int atNum2)
   {
      if( atNum2 != BROMINE )
         return at1.GetHvyDegree() == at2.GetHvyDegree() ? 0.2D : 0D;

      /////// BROMINE BROMINE
      if( at1.GetValence() == at2.GetValence() )
         return 1D;

      return 0.0D;
   }
}
