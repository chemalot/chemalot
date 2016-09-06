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

package com.gNova.circularFP;

import java.util.*;

public interface Fingerprint
{
  /**
    * hashed integer atom ID's
    */
   public String getAtomIDString();

  /**
    * bit0 is on far left.
    */
   public String getBitString();

   /**
    * bit0 is on far left.
    *
    * This will always have an even number of characters (last will be filled
    * with 0 if needed.
    */
   public String getHexString();

   /**
    * index of bits set in fingerprint ordered ascending.
    */
   public int[] getBits();

   public double Tanimoto(Fingerprint other);

   /**
    * Return the number of the bits set.
    */
   public int getNBits();

   /**
    * Return the Bit Set
    */
   public BitSet getBitSet();

   /**
    * Create a folded fingerprint of the specified number of bits by using modulus.
    * @return a new folded fingerprint instance.
    */
   //public Fingerprint fold(int size);
}
