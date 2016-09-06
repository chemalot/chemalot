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

public class mean extends AggFunction
{  public static final String DESCRIPTION
         = "[outName[:format] = ] mean([distinct] Tag): returns the mean of tag 'Tag'.\n";



   public mean(String outTag, String resultFormat, String funcName, String funcArg)
   {  super(outTag, resultFormat == null ? "si2" : resultFormat, funcName, funcArg);
   }



   private String getResult(String fmtStr)
   {  if( valueContainer.size() == 0 )
         return "";

      double sum = 0D;
      for(String v : valueContainer)
      {  try
         {  sum += Float.parseFloat(v);
         } catch (NumberFormatException e)
         {  System.err.printf("Warning: %s for %s cannot be converted to float.\n",
                              v, getOutTagName());
         }
      }

      return DataFormat.formatNumber(sum/valueContainer.size(), fmtStr);
   }

   @Override
   public String getResult(int indxInGrp)
   {  return getResult(resultFormat);
   }
}
