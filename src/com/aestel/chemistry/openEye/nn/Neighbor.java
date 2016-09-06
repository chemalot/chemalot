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
package com.aestel.chemistry.openEye.nn;

/**
 * For sorting neighbors by descendign similarity
 * @author albertgo
 *
 */
class Neighbor implements Comparable<Neighbor>
{  final int neighBorIdx;
   final double neighBorSim;

   public Neighbor(int i, double sim)
   {  neighBorIdx = i;
      neighBorSim = sim;
   }

   @Override
   public int compareTo(Neighbor o)
   {  int v = Double.compare(o.neighBorSim, neighBorSim);
      if( v == 0 )
         return neighBorIdx - o.neighBorIdx;
      return v;
   }

   public boolean equals(Neighbor o)
   {  return neighBorIdx == o.neighBorIdx && neighBorSim == o.neighBorSim;
   }
}
