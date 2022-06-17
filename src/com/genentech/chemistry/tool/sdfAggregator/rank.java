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
import java.util.Iterator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.aestel.utility.ArrayIndexComparator;

/**
 * This class implements the count function.
 * @author Johnny Wu Aug 11, 2011
 *
 */


class rank extends AggFunction
{  public static final String DESCRIPTION
      = "[outName = ] rank([distinct] [asc|desc] fieldName): returns the rank of the value in the group use with outputAll.\n";

   private static final Pattern ASCDesc = Pattern.compile(
                        "\\s*((asc|desc)\\s+)?(.+)", Pattern.CASE_INSENSITIVE);
   private final double multiplier;

   public rank(String outTag, String resultFormat, String funcName, String funcArg)
   {  super(outTag, resultFormat, funcName, funcArg);

      if( resultFormat != null ) throw new Error("Formattign not supported for concatenate");

      Matcher mat = ASCDesc.matcher(funcArg);
      if( !mat.matches() )
         throw new Error(String.format("Invalid argument to rank: %s", funcArg));
      multiplier = "desc".equalsIgnoreCase(mat.group(2)) ? -1D : 1D;
   }

   @Override
   protected String getTagName(String funcArg)
   {  Matcher mat = ASCDesc.matcher(funcArg);
      if( !mat.matches() )
         throw new Error(String.format("Invalid argument to rank: %s", funcArg));

      return mat.group(3);
   }

   @Override
   public String getResult(int indxInGrp)
   {  if( valueContainer.size() == 0 ) return "";
   
      int[] ranks = getRanks();

      return Integer.toString(ranks[indxInGrp]);
   }

   
   protected int[] getRanks()
   {  Double[] valD = new Double[valueContainer.size()];
      Iterator<String> it = valueContainer.iterator();

      for( int i=0; it.hasNext(); i++)
      {  String val = it.next();
         if(val.length() > 0)
            valD[i] = Double.valueOf(val) * multiplier;
         else if( multiplier == 1D )
            valD[i] = Double.NEGATIVE_INFINITY;
         else
            valD[i] = Double.POSITIVE_INFINITY;
      }

      ArrayIndexComparator<Double> aic = new ArrayIndexComparator<Double>(valD);
      Integer[] indexes = aic.createIndexArray();
      Arrays.sort(indexes, aic);

      // compute rank
      int[] ranks = new int[indexes.length];
      int rank = 0 ;
      double lastEle = indexes.length > 0 ? valD[indexes[0]] : 0D;
      for(int i=0; i< indexes.length; i++)
      {  if( valD[indexes[i]] != lastEle )
            rank++;
         ranks[indexes[i]] = rank;
      }
      return ranks;
   }
}
