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
package com.genentech.retrival.SDFExport;

import java.io.*;

public class SDFParser
{

   private BufferedReader in;
   private StringBuilder mol;
   private StringBuilder data;
   private int nStruct = 0;
   private SDRecord lastSDRec;

   public SDFParser(File sdFile) throws FileNotFoundException
   {  if (sdFile == null)
         in = new BufferedReader(new InputStreamReader(System.in));
      else
         in = new BufferedReader(new FileReader(sdFile));
      mol = new StringBuilder();
      data = new StringBuilder();
   }

   public boolean hasNext() throws IOException
   {  if (lastSDRec != null)
         return true;

      mol.setLength(0);
      data.setLength(0);
      String line;
      boolean inMolFile = true;

      while (true)
      {
         if ((line = in.readLine()) == null)
         {  if ("".equals(mol.toString().trim()))
               return false;

            if (!inMolFile)
               throw new Error("Invalid end of sd file!");

            line = "$$$$";
         }

         if (line.startsWith("$$$$"))
         {  nStruct++;
            lastSDRec = new SDRecord(mol.toString(), data.toString());
            return true;
         } else if (!inMolFile || line.startsWith(">"))
         {  inMolFile = false;
            data.append(line).append("\n");
         } else
         {  mol.append(line).append("\n");
         }
      }
   }

   public SDRecord next()
   {  SDRecord sdRec = lastSDRec;
      lastSDRec = null;
      return sdRec;
   }
}
