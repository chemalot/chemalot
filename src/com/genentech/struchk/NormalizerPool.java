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

package com.genentech.struchk;

import java.util.HashSet;
import java.util.Set;

import openeye.oechem.OEGraphMol;

import org.apache.commons.pool.BasePoolableObjectFactory;
import org.apache.commons.pool.impl.GenericObjectPool;
import org.apache.commons.pool.impl.GenericObjectPool.Config;

import com.aestel.utility.LogHelper;
import com.genentech.oechem.tools.OETools;
import com.genentech.struchk.oeStruchk.StruChkHelper.CHECKType;

/**
 * Pool of {@link Normalizers} objects.
 *
 * @author A. Gobbi 2006 Copyright 2006 Genentech Inc.
 */
public class NormalizerPool
{
   public static final NormalizerPool DEFAULT_POOL = new NormalizerPool(
         new HashSet<CHECKType>(0), false);

   private final GenericObjectPool pool;

   private NormalizerPool(Set<CHECKType> exclusions, boolean errorAsWarning)
   {
      Config poolcfg = new Config();
      poolcfg.maxIdle = 3;
      poolcfg.maxActive = 20;
      poolcfg.timeBetweenEvictionRunsMillis = 1000 * 3600;
      poolcfg.testWhileIdle = true;
      // since we do not set maxActive for the moment this does not matter
      poolcfg.whenExhaustedAction = GenericObjectPool.WHEN_EXHAUSTED_BLOCK;
      poolcfg.maxWait = 40000;
      NormalizerFactory fact = new NormalizerFactory(exclusions, errorAsWarning);

      pool = new GenericObjectPool(fact, poolcfg);
   }

   /**
    * Check a molfile.
    * @param gneStructFlag may be "" in which case the structure flag is guessed from the
    *        specified stereo centers.
    */
   public GNEMolecule normalizeMol(String molStr, String gneStructFlag, int nNonTetrahedralChiral)
   {
      PooledNormalizer pnorm = null;
      try
      {
         pnorm = (PooledNormalizer) pool.borrowObject();
         return pnorm.norm.normalizeMol(molStr, gneStructFlag, nNonTetrahedralChiral);
      } catch (Exception e)
      {
         handleException(pnorm, e);
         throw new Error(e);
      } finally
      {
         returnToPool(pnorm);
      }
   }


   public String validateMol(String molStr, String structFlag, int nNonTetrahedralChiral,
         String stereoComment)
   {

      GNEMolecule gneMol = normalizeMol(molStr, structFlag, nNonTetrahedralChiral);
      if (gneMol.hasError())
         return gneMol.getErrors();

      if (!structFlag.equals("No Stereo")
            && !structFlag.equals("Single Known Stereoisomer")
            && !(gneMol.getNChiral() == 1 && structFlag
                  .equals("Mixture of Enantiomers"))
            && (stereoComment == null || stereoComment.length() == 0))
         return "Please provide a structure comment with information about the unspecified stereo centers.\n";

      return null;
   }

   public GNEMolecule normalizeSmi(String smi, String gneStructFlag, int nNonTetrahedralChiral)
   {
      PooledNormalizer pnorm = null;
      try
      {
         pnorm = (PooledNormalizer) pool.borrowObject();
         return pnorm.norm.normalizeSmi(smi, gneStructFlag, nNonTetrahedralChiral);
      } catch (Exception e)
      {
         handleException(pnorm, e);
         throw new Error(e);
      } finally
      {
         returnToPool(pnorm);
      }
   }

   /**
    * Check a molfile.
    */
   public GNEMolecule normalizeOEMol(OEGraphMol mol, String gneStructFlag, int nNonTetrahedralChiral)
   {
      PooledNormalizer pnorm = null;
      try
      {
         pnorm = (PooledNormalizer) pool.borrowObject();
         return pnorm.norm.normalizeOEMol(mol, gneStructFlag, nNonTetrahedralChiral);
      } catch (Exception e)
      {
         handleException(pnorm, e);
         throw new Error(e);
      } finally
      {
         returnToPool(pnorm);
      }
   }

   /**
    * Assumes 0 nNonTetrahedralChiral
    *
    * @deprecated use {@link #normalizeMol(String, String, int)}
    */
   @Deprecated
   public GNEMolecule normalizeMol(String molStr, String gneStructFlag)
   {  return normalizeMol(molStr, gneStructFlag, 0);
   }


   /**
    * Assumes 0 nNonTetrahedralChiral
    *
    * @deprecated use {@link #validateMol(String, String, int, String)}
    */
   @Deprecated
   public String validateMol(String molStr, String structFlag,
         String stereoComment)
   {  return validateMol(molStr, structFlag, 0, stereoComment);
   }


   /**
    * Assumes 0 nNonTetrahedralChiral
    *
    * @deprecated use {@link #normalizeSmi(String, String, int)}
    */
   @Deprecated
   public GNEMolecule normalizeSmi(String smi, String gneStructFlag)
   {  return normalizeSmi(smi, gneStructFlag, 0);
   }

   /**
    * Assumes 0 nNonTetrahedralChiral
    *
    * @deprecated use {@link #normalizeSmi(String, String, int)}
    */
   @Deprecated
   public GNEMolecule normalizeOEMol(OEGraphMol mol, String gneStructFlag)
   {  return normalizeOEMol(mol, gneStructFlag, 0);
   }



   public void returnToPool(PooledNormalizer norm)
   {
      try
      {
         pool.returnObject(norm);
      } catch (Exception e)
      {
         LogHelper.severe(e); // this is finally block there is not much more we
                              // can do
      }
   }

   private static void handleException(PooledNormalizer pnorm, Exception e)
   {
      LogHelper.severe(e);
      if (pnorm != null)
      { // nothing to do for now
      }
   }

   static class PooledNormalizer
   {
      public PooledNormalizer(Normalizer norm)
      {  assert norm != null;
         this.norm = norm;
      }

      final Normalizer norm;
      final long checkOutTime = System.currentTimeMillis();
   }

   static class NormalizerFactory extends BasePoolableObjectFactory
   {
      private final boolean errorAsWarning;
      private final Set<CHECKType> exclusions;

      public NormalizerFactory(Set<CHECKType> exclusions, boolean errorAsWarning)
      {
         this.exclusions = exclusions;
         this.errorAsWarning = errorAsWarning;
      }

      /**
       * Creates one instance of SmiTalk for the pool.
       */
      @Override
      public Object makeObject() throws Exception
      {
         return new PooledNormalizer(new Normalizer(exclusions, errorAsWarning));
      }

      @Override
      public void passivateObject(Object obj) throws Exception
      { // not used
      }

      @Override
      public void destroyObject(Object obj)
      {
         ((PooledNormalizer) obj).norm.close();
      }

      @Override
      public boolean validateObject(Object obj)
      {
         if (System.currentTimeMillis() - ((PooledNormalizer) obj).checkOutTime < 1000 * 3600 * 20)
            return false; // older than 12 hours
         return true;
      }
   }

   public static void main(String... args)
   {
      String smi = "C[C@H](F)O";
      OEGraphMol mol = new OEGraphMol();
      OETools.smiToMol(mol, smi);
      String molStr = OETools.molToString(mol);
      mol.delete();

      String stereo = "Single Unknown Stereoisomer";
      stereo = "";
      System.err.println(NormalizerPool.DEFAULT_POOL.validateMol(molStr,
            stereo, "hello"));
      System.err.println(NormalizerPool.DEFAULT_POOL.normalizeMol(molStr,
            stereo).getTautomerISmi());

   }
}
