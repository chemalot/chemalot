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

import com.aestel.utility.DataFormat;

public class stddev extends CachedAggFunction
{  public static final String DESCRIPTION
         = "[outName[:format] = ] stddev([distinct] Tag): returns the stddev of tag 'Tag'.\n";



   public stddev(String outTag, String resultFormat, String funcName, String funcArg)
   {  
      super(outTag, resultFormat == null ? "si2" : resultFormat, funcName, funcArg);
   }
   
   
   @Override
   String getResult(String fmtStr)
   {  
      if( valueContainer.size() < 2 )
         return "";

      double sPrev;
      double s = 0D;
      double var = 0D;
      double avgPrev;
      double avg =0D;
      int cntr = 0;
      for(String v : valueContainer)
      {  try
         {  double vd = Double.parseDouble(v);
            ++cntr;

           if (cntr == 1)
           {   // Set the very first values.
               avg = vd;
           } else
           {   // Save the previous values.
               avgPrev = avg;
               sPrev   = s;

               // Update the current values.
               avg = avgPrev + (vd - avgPrev) / cntr;
               s   = sPrev   + (vd - avgPrev) * (vd - avg);
               var = s / (cntr - 1);
           }
         } catch (NumberFormatException e)
         {  System.err.printf("Warning: %s for %s cannot be converted to float.\n",
                              v, getOutTagName());
         }
      }

      // Return sqrt of variance = stdev
      return DataFormat.formatNumber(Math.sqrt(var), fmtStr);
   }
}
