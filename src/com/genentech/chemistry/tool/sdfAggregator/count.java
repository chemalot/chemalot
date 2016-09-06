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
 * This class implements the count function.
 * @author Johnny Wu Aug 11, 2011
 *
 */


class count extends AggFunction
{  public static final String DESCRIPTION
      = "[outName = ] count([distinct] fieldName): returns the number records having a value in fieldName.\n";

   public count(String outTag, String resultFormat, String funcName, String funcArg)
   {  super(outTag, resultFormat, funcName, funcArg);
      if( resultFormat != null ) throw new Error("Formattign not supported for concatenate");
   }


   @Override
   public String getResult (int indxInGrp)
   {  return Integer.toString(valueContainer.size());
   }

}
