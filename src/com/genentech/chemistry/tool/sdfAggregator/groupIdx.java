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
package com.genentech.chemistry.tool.sdfAggregator;

/**
 * This class implements the max function.
 * @author A Gobbi 2013
 *
 */

import openeye.oechem.OEGraphMol;

public class groupIdx  implements AggInterface
{  public static final String DESCRIPTION
         = "outName = groupIdx(): returns index of each group defiend by the groupby in order of input.\n";

   private final String outTag;

   private int counter = 0;

   public groupIdx(String outTag, String resultFormat, String funcName, String funcArg)
   {  this.outTag  = outTag.trim();
      assert resultFormat == null || resultFormat.trim().length() == 0 : "no foramt supported for counter";
      assert funcArg == null || funcArg.trim().length() == 0 : "No column supported for counter";
   }


   @Override
   public String getResult(int indxInGrp)
   {  return Integer.toString(counter);
   }


   @Override
   public void init()
   {  counter++;
   }


   @Override
   public void process(OEGraphMol mol)
   {  //nothing to do
   }


   @Override
   public String getOutTagName()
   {  return outTag;
   }
}
