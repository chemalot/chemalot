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
 * This interface provides the definition for aggregation functions.
 * @author Johnny Wu Aug 11, 2011
 *
 */

import openeye.oechem.OEGraphMol;

interface AggInterface {
   public void init ();
   public void process (OEGraphMol mol);
   public String getOutTagName ();
   public String getResult(int indxInGrp);

   /**
    * In addition an aggregation function must implement a constructor with a
    * Single String argument which is called containing the arguments passed.
    *
    *  The constructor may throw an IllegalArgumentException
    *
    *  public AggInterface(String outTag, String funcName, String funcArg) throws IllegalArgumentException
    */


   /**
    * In addition an aggregation function must implement:
    * static public final String DESCRIPTION
    */

}
