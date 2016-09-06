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
 * This class implements the standard deviation function.
 * @author Ben Sellers Feb, 2015
 *
 */

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import com.aestel.utility.DataFormat;

public class stddev extends AggFunction
{  public static final String DESCRIPTION
         = "[outName[:format] = ] stddev([distinct] Tag): returns the stddev of tag 'Tag'.\n";



   public stddev(String outTag, String resultFormat, String funcName, String funcArg)
   {  
      super(outTag, resultFormat == null ? "si2" : resultFormat, funcName, funcArg);
   }
   
   private String getResult(String fmtStr)
   {  
      if( valueContainer.size() == 0 )
         return "";

      List<Double> numericValues = new ArrayList<Double>();
      Iterator<String> it = valueContainer.iterator();

      //
      // Parse
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
      
      //
      // Calc mean
      double sum = 0D;
      for(Double value : numericValues)
      {     
         sum += value;
      }
      double mean = sum/numericValues.size(); 
      
      //
      // Calculate the variance
      double sumOfSquaredDiffs = 0D;
      double variance          = 0D;      
      for(Double value : numericValues)
      {       
         sumOfSquaredDiffs += (value - mean) * (value - mean);
      }
      
      variance = sumOfSquaredDiffs / numericValues.size();

      // Return sqrt of variance = mstdev
      return DataFormat.formatNumber(Math.sqrt(variance), fmtStr);
   }
   

   @Override
   public String getResult(int indxInGrp)
   {  return getResult(resultFormat);
   }
}
