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

/**
 * Factory class to create comparator and comparable objects.
 *
 * The comparator can calculate its similarity to a comparable. But not vice versa.
 *
 * @author albertgo
 *
 * @param <IN> the input object which is can be converted into a comparable.
 * @param <T> the comparable class type
 * @param <Y> the comparator class type
 */
public interface SimComparatorFactory<IN, T, Y extends SimComparator<T>>
{  Y createComparator(T comparable);

   T createComparable(IN in);

   void close();
}
