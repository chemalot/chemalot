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

import java.util.ArrayList;
import java.util.List;

import openeye.oechem.OEIFlavor;
import openeye.oechem.OEMDLQueryOpts;
import openeye.oechem.OEMolBase;
import openeye.oechem.OEQMol;
import openeye.oechem.OESubSearch;
import openeye.oechem.oechem;
import openeye.oechem.oemolistream;

public class MDLSSSMatcher
{  private final List<OESubSearch> queries;
   private final List<String>      qIds;

   /**
    * @param qFile opened with OEAroModelOpenEye, so default streams should work fine.
    */
   public MDLSSSMatcher(String qFile)
   {  queries = new ArrayList<OESubSearch>(200);
      qIds    = new ArrayList<String>(200);

      // read queries into List
      oemolistream qfile = new oemolistream(qFile);
      int aromodel = OEIFlavor.Generic.OEAroModelOpenEye;
      int qflavor  = qfile.GetFlavor(qfile.GetFormat());
      qfile.SetFlavor(qfile.GetFormat(),(qflavor|aromodel));
      int opts = OEMDLQueryOpts.Default|OEMDLQueryOpts.SuppressExplicitH;
      OEQMol qmol = new OEQMol();
      int nMol = 0;
      while( oechem.OEReadMDLQueryFile(qfile,qmol,opts) )
      {  OESubSearch ss = new OESubSearch(qmol);

         queries.add(ss);

         String qId = qmol.GetTitle(); // note: oechem.OEGetSDData() does not work on qMols
         if( qId == null ||qId.length() == 0)
         {  System.err.printf("Molecule %d has no title, nondeterministic id %d will be used!\n",
                     nMol, nMol);
            qId = Integer.toString(nMol);
         }

         qIds.add(qId);
         nMol++;
         qmol.Clear();
      }
      System.err.printf("Read %d queriy molecules\n", nMol);

      qmol.delete();
      qfile.close();
      qfile.delete();
   }

   /**
    * Calling this method is threads safe, no internal state is changed.
    */
   public boolean findMatches(OEMolBase mol, boolean firstMatch)
   {  StringBuilder sb = new StringBuilder(50);

      for(int i=0; i<queries.size(); i++)
      {  OESubSearch q = queries.get(i);

         if(q.SingleMatch(mol))
         {  sb.append(qIds.get(i));
            if( firstMatch ) break;

            sb.append(".");
         }
      }
      if(sb.length() > 0)
      {  if( ! firstMatch ) sb.setLength(sb.length()-1); // remove trailing "."

         oechem.OESetSDData(mol, "SSSMatchName", sb.toString());
         return true;
      }
      return false;
   }

   public void close()
   {  for( OESubSearch q : queries )
         q.delete();
   }
}
