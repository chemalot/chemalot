/*
   Copyright 2006-2014 Man-Ling Lee & Alberto Gobbi

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

Contact: aestelSW@gmail.com
*/

package com.aestel.chemistry.openEye;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

import openeye.oechem.OEGraphMol;
import openeye.oechem.OEMolBase;
import openeye.oechem.oechem;
import openeye.oechem.oemolithread;
import openeye.oechem.oemolothread;

/**
 * SphereExclusion method implementation with generic comparable and comparator types.
 *
 * @author albertgo
 *
 */
public class SphereExclusion<T, Y extends SimComparator<T>>
{  private final oemolothread ofs;
   private final double radius;
   private final boolean reverseMatch;
   private final boolean printAll;
   private final boolean printSphereMatchCount;
   private final List<Y> centroids;
   private final SimComparatorFactory<OEMolBase, T, Y> comparableFact;
   private long start;
   private int iCounter;
   private int exCounter;
   private int incCounter;

   /**
    *
    * @param comparableFact factory class to create the comparator.
    * @param refFile filename of file containing pre-existing reference compounds (centroids)
    * @param outFile output filename
    * @param radius  molecules with sim >= radius are considered cluster members
    * @param reverseMatch compare candidates to centroids starting with the latest
    *        found centroid, this is faster than using the first centroid, when
    *        using printSphereMatchCount reverseMatch false is more useful.
    * @param printSphereMatchCount check and print membership of candidates not
    *        only to the first centroid which has sim >= radius but to all centroids
    *        found up to that input. This will output a candidate multiple times.
    * @param printAll print also members not only centroids
    */
   public SphereExclusion(SimComparatorFactory<OEMolBase, T, Y> comparableFact,
                    String refFile, String outFile, double radius,
                    boolean reverseMatch, boolean printSphereMatchCount, boolean printAll)
   {  this.comparableFact = comparableFact;
      this.ofs = new oemolothread(outFile);
      this.radius = radius;
      this.reverseMatch = reverseMatch;
      this.printAll = printAll;
      this.printSphereMatchCount = printSphereMatchCount;
      this.centroids = new ArrayList<Y>(2000);
      if( refFile != null ) readReferenceFile(refFile);
   }


   public void run(String inFile)
   {  oemolithread ifs = new oemolithread(inFile);

      start = System.currentTimeMillis();
      iCounter = 0;
      exCounter = 0;
      incCounter = 0;

      OEMolBase mol = new OEGraphMol();
      while(oechem.OEReadMolecule(ifs, mol))
      {  iCounter++;

         T tmp = comparableFact.createComparable(mol);
         Y current = comparableFact.createComparator(tmp);
         // tmp needs delete;

         double maxSim = -1D;
         int sphereMatchCounter=0;
         int centIdx;
         final int increment;
         final int lastCentroid;

         if( reverseMatch )
         {  // compare next record with last centroid found first
            // this is usually faster since it is assumed that close input
            // records are more similar to each other and the likelihood of
            // Exclusion is therefore higher.
            increment = -1;
            lastCentroid = -1;
            centIdx = centroids.size();
         }else
         {  // compare next record with first centroid found first
            // this assigns the new record to the first matching centroid which
            // might have superior properties.
            increment = 1;
            lastCentroid = centroids.size();
            centIdx = -1;
         }

         // compare to all so far known centroids
         while( (centIdx += increment) != lastCentroid )
         {  Y c2 = centroids.get(centIdx);
            double sim = c2.similarity(current);
            if( sim > maxSim ) maxSim = sim;
            if( sim >= radius )  // inside the radius of centrIdx's sphere?
            {  exCounter++;
               sphereMatchCounter++;
               oechem.OESetSDData(mol, "sphereIdx", Integer.toString(centIdx));
               oechem.OESetSDData(mol, "maxSim", String.format("%.3f",maxSim));

               if( printSphereMatchCount )
               {  oechem.OESetSDData(mol, "sphereMatchCounter", Integer.toString(sphereMatchCounter));
                  oechem.OEWriteMolecule(ofs, mol);

               } else // do not assign to multiple clusters
               {  if( printAll ) oechem.OEWriteMolecule(ofs, mol);
                  break;
               }
            }
         }

         if( sphereMatchCounter == 0 ) // was not a member of any centroids
         {  // so it becomes a new centroid
            String spherIdx = Integer.toString(centroids.size());
            oechem.OESetSDData(mol, "sphereIdx",  spherIdx);
            oechem.OESetSDData(mol, "includeIdx", spherIdx);
            oechem.OESetSDData(mol, "maxSim",     "1");
            centroids.add(current);
            incCounter++;

            oechem.OEWriteMolecule(ofs, mol);
         }

         if(iCounter % 100 == 0) System.err.print(".");
         if(iCounter % 4000 == 0)
         {  System.err.printf( " %d %d included %dsec\n",
                  iCounter, incCounter, (System.currentTimeMillis()-start)/1000);
         }
      }

      ifs.close();
      inFile = inFile.replaceAll(".*" + Pattern.quote(File.separator), "");
      System.err.printf("SphereExclusion: Read %d structures from %s, %d included, %d excluded. %d sec\n",
            iCounter, inFile, incCounter, exCounter, (System.currentTimeMillis()-start)/1000);
   }

   /**
    * Read in file with preselected centroids
    */
   private void readReferenceFile(String refFile)
   {  OEMolBase mol = new OEGraphMol();
      oemolithread ifs = new oemolithread(refFile);

      while(oechem.OEReadMolecule(ifs, mol))
      {  iCounter++;
         Y c = comparableFact.createComparator(comparableFact.createComparable(mol));
         centroids.add(c);
      }
      System.err.printf("%d reference molecules read.\n", iCounter);
      ifs.close();
      mol.delete();
      ifs.delete();
   }


   public void close()
   {  ofs.close();
      ofs.delete();
      for(SimComparator<?> c: centroids)
         c.close();

      comparableFact.close();
   }
}

