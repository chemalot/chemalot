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

import java.util.Arrays;
import java.util.Collections;
import java.util.stream.Collectors;

import com.aestel.utility.DataFormat;

/**
 * This class implements the rankPercent function.
 * @author Alberto Gobbi 8/2021
 *
 */


class rankPct extends rank
{  @SuppressWarnings("hiding")
   public static final String DESCRIPTION
      = "[outName = ] rankPct([distinct] [asc|desc] fieldName): returns the percentage rank of the value in the group (use with outputAll).\n";

   public rankPct(String outTag, String resultFormat, String funcName, String funcArg)
   {  super(outTag, resultFormat, funcName, funcArg);
   }

   @Override
   public String getResult(int indxInGrp)
   {  if( valueContainer.size() == 0 ) return "";
      if( valueContainer.size() == 1 ) return "0.5";
   
      int[] ranks = getRanks();
      int max = Collections.max(Arrays.stream(ranks).boxed().collect(Collectors.toList()));
      if( max == 0 ) max++;
      
      return DataFormat.formatNumber(ranks[indxInGrp]/(double)(max) * 100.,"r3");
   }
}
