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

/**
 * A Fingerprint encodes the presence or absence of features.
 *
 * @author albertgo
 *
 */
public interface Fingerprint
{  /**
    * bit0 is on far left.
    * @return string of 0 and 1.
    */
   public String getBinString();

   /**

    * bit0 is on far left.
    *
    * This will always have an even number of characters (last will be filled
    * with 0 if needed.
    *
    * @return hexadecimal encoding of fingeprint
    */
   public String getHexString();

   /**
    * @return index of bits set to 1 in fingerprint ordered ascending.
    */
   public int[] getBits();

   /**
    * Compute similarity using tanimoto coefficient.
    */
   public double tanimoto(Fingerprint other);
//   public double Tversky(Fingerprint other, double a, double b);

   /**
    * Return the number of the bits set to 1.
    */
   public int getNBits();

   /**
    * Create a folded fingerprint of the specified number of bits by using modulus.
    * @return a new folded fingerprint instance.
    */
   public Fingerprint fold(int size);
}
