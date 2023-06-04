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
 * This class implements the mean function.
 * @author Johnny Wu Aug 11, 2011, AG 2012
 *
 */

import com.aestel.utility.DataFormat;

public class mean extends CachedAggFunction
{  public static final String DESCRIPTION
         = "[outName[:format] = ] mean([distinct] Tag): returns the mean of tag 'Tag'.\n";


   public mean(String outTag, String resultFormat, String funcName, String funcArg)
   {  super(outTag, resultFormat == null ? "si2" : resultFormat, funcName, funcArg);
   }


   @Override
   String getResult(String fmtStr)
   {  if( valueContainer.size() == 0 )
         return "";

      double avg = 0D;
      int cntr = 0;
      for(String v : valueContainer)
      {  try
         {  double vd = Double.parseDouble(v);
            avg += (vd-avg)/++cntr;
         } catch (NumberFormatException e)
         {  System.err.printf("Warning: %s for %s cannot be converted to float.\n",
                              v, getOutTagName());
         }
      }

      return DataFormat.formatNumber(avg, fmtStr);
   }
}
