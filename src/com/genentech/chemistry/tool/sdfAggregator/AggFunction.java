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
import java.util.Collection;
import java.util.HashSet;

import openeye.oechem.OEGraphMol;
import openeye.oechem.oechem;

/**
 * This abstract class provides the definition of aggregation functions.
 * @author Johnny Wu Aug 11, 2011, AG 2012
 *
 */
abstract class AggFunction implements AggInterface
{  protected final boolean distinct;
   protected final String tagName;
   protected final Collection<String> valueContainer;
   protected final String outTag;
   protected final String resultFormat;

   public AggFunction(String outTag, String resultFormat, String funcName, String funcArg)
   {  this.outTag  = outTag.trim();
      this.resultFormat = resultFormat == null ? null : resultFormat.trim();

      if( "distinct".equals(funcArg.trim() ) )
      {  distinct = true;
         tagName = null;
         valueContainer = new HashSet<String>();
         return;
      }

      if( funcArg.trim().startsWith("distinct "))
      {  distinct = true;
         funcArg = funcArg.substring(funcArg.indexOf("distinct ")+9);
      } else
      {  distinct = false;
      }

      tagName = getTagName(funcArg);

      if( distinct )
         valueContainer = new HashSet<String>();
      else
         valueContainer = new ArrayList<String>();
   }

   /** to be overwritten to parse more complicated parameters */
   protected String getTagName(String funcArg)
   {  return funcArg.trim();
   }


   @Override
   public void process (OEGraphMol mol)
   {  String tagVal = oechem.OEGetSDData(mol, tagName);
      if( tagVal != null && tagVal.length() > 0)
         valueContainer.add(tagVal);
   }

   @Override
   public void init()
   {  valueContainer.clear();
   }

   public String getAggregatedFieldName()
   {  return tagName;
   }

   @Override
   public String getOutTagName ()
   {  return outTag;
   }


}
