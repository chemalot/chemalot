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
 * This class implements the boltzman probability computation form relative energy.
 * @author AG 11, 2021
 *
 */

import com.aestel.utility.DataFormat;

public class boltzmannProbability extends AggFunction
{  public static final String DESCRIPTION
         = "[outName[:format] = ] boltzmannProbability(relEnergyKCalMol_Tag): returns the probability of this record given the energy relative to other energies (T=298)\n";

   protected final double RT = 0.5921856D;  // at 298K


   public boltzmannProbability(String outTag, String resultFormat, String funcName, String funcArg)
   {  super(outTag, resultFormat == null ? "si2" : resultFormat, funcName, funcArg);
   }


   @Override
   public String getResult(int indxInGrp)
   {  if( valueContainer.size() == 0 )
         return "";
   
      double partitionFunc = 0;
      
      // Get partition function
      for(String v : valueContainer)
      {  try
         {  double deltaE = Double.parseDouble(v);
            double bFactor = Math.exp(-1 * deltaE / RT);
            partitionFunc+=bFactor;
         } catch (NumberFormatException e)
         {  System.err.printf("Warning: %s for %s cannot be converted to float.\n",
                              v, getOutTagName());
         }
      }
   
      // Get current energy
      double bFactor = 0;
      try
      {  double deltaE = Double.parseDouble((String)valueContainer.toArray()[indxInGrp]);
         bFactor = Math.exp(-1 * deltaE / RT);
      } 
      catch (NumberFormatException e)
      {  // Msg already printed above
         return "";
      }
      
   
      return DataFormat.formatNumber(bFactor/partitionFunc, resultFormat);
   }
}
