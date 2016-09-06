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
package com.aestel.chemistry.openEye;

public interface SimComparator<T>
{  double similarity(T other);

   /**
    * Compare this T with the other and return the similarity.
    *
    * This method needs only to support comparisons with the specific implementation.
    * eg. {@link com.aestel.chemistry.openEye.fp.FPComparator} only needs to implement similarity(FPComparator other).
    *
    * A {@link ClassCastException} may be thrown if a subtype is passed.
    *
    */
   double similarity(SimComparator<T> other);

   /**
    * Compare this T with the other and return the similarity.
    *
    * This method is allowed to return 0 if the similarity is below minSim
    *
    */
   double similarity(SimComparator<T> other, double minSim);

   /**
    * This should delete any OE objects cloned in the constructor.
    */
   void close();
}
