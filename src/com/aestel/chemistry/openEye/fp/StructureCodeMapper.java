/*
   Copyright 2006-2014 Man-Ling Lee & Alberto Gobbi

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

Contact: aestelSW@gmail.com
*/

package com.aestel.chemistry.openEye.fp;

import java.util.Iterator;

/**
 * Map a structural feature to an bit position in a fingerprint.
 * @author albertgo
 *
 */
public interface StructureCodeMapper
{  public static final int TYPEIdx = 0;
   public static final int NAMEIdx = 1;
   public static final int INDEXIdx = 2;

   /** @return bit position of feature fragName. */
   public int getIndex(String type, String fargName);

   /** @return name of feature with position idx */
   public String getName(int idx);

   public int getMinIdx();
   public int getMaxIdx();

   public void close();

   /**
     * @return iterator of String[3] = {CodeType, CodeName, BitIndex}.
    */
   public Iterator<String[]> getIterator();
}
