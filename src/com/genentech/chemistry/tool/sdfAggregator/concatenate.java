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


class concatenate extends AggFunction
{  public static final String DESCRIPTION
      = "[outName = ] concatenate([distinct] fieldName[ , separator]): returns the concatenated list of values in fieldName.\n";

   private final String sep;

   public concatenate(String outTag, String resultFormat, String funcName, String funcArg)
   {  super(outTag, resultFormat, funcName, funcArg);
      if( resultFormat != null ) throw new Error("Formattign not supported for concatenate");
      if( funcArg.contains(",") )
         sep = funcArg.substring(funcArg.indexOf(',')+1);
      else
         sep = ", ";
   }

   @Override
   protected String getTagName(String funcArg)
   {  if( funcArg.contains(",") )
         funcArg = funcArg.substring(0,funcArg.indexOf(','));
      return super.getTagName(funcArg);
   }

   @Override
   public String getResult(int indxInGrp)
   {  StringBuilder sb  = new StringBuilder(valueContainer.size()*10);
      for( String val: valueContainer)
         sb.append(val).append(sep);

      if(sb.length()>0) sb.setLength(sb.length()-sep.length());

      return sb.toString();
   }

}
