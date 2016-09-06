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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;




import java.util.List;

import com.aestel.utility.ArrayIndexComparator;
import com.aestel.utility.DataFormat;

/**
 * This class implements the median function.
 * @author Ben Sellers Feb. 2015
 */


class median extends AggFunction
{  public static final String DESCRIPTION
      = "[outName[:format] = ] median([distinct] tagName): returns the median of tagName values.\n";

   public median(String outTag, String resultFormat, String funcName, String funcArg)
   {  
      super(outTag, resultFormat == null ? "si2" : resultFormat, funcName, funcArg);
   }


   private String getResult(String fmtStr)
   {  
      if( valueContainer.size() == 0 )
         return "";
      
      double median = 0D;
      List<Double> numericValues = new ArrayList<Double>();
      Iterator<String> it = valueContainer.iterator();

      String val = "";
      try
      {
         while (it.hasNext())
         {  
            val = it.next();
            if(val.trim().length() > 0)
               numericValues.add(Double.valueOf(val.trim()));
         }
      }
      catch (NumberFormatException e)
      {  
         System.err.printf("Warning: %s for %s cannot be converted to float.\n",
                  val, getOutTagName());
      }
   
      Double [] numericValuesArray = numericValues.toArray(new Double [numericValues.size()]);
      ArrayIndexComparator<Double> aic = new ArrayIndexComparator<Double>(numericValuesArray);
      Integer[] indexes = aic.createIndexArray();
      Arrays.sort(indexes, aic);

      // compute median
      if (indexes.length % 2 == 0)
      {
         // even -- take mean of two middle values
         int nLowerMiddleIndex = ((indexes.length)/2) - 1;
         int nUpperMiddleIndex = nLowerMiddleIndex+1;
         median = (numericValues.get(indexes[nLowerMiddleIndex]) + numericValues.get(indexes[nUpperMiddleIndex]) ) / 2;
      }
      else
      {
         // odd -- take middle value
         int nMiddleIndex = ((indexes.length+1)/2) - 1;
         median = numericValues.get(indexes[nMiddleIndex]);
      }

      return DataFormat.formatNumber(median, fmtStr);
   }
   
   @Override
   public String getResult(int indxInGrp)
   {  return getResult(resultFormat);
   }
}
