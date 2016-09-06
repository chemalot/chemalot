/*
   Copyright 2008-2014 Genentech Inc.

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
package com.genentech.chemistry.openEye;

/** computes the similarity of the head atoms for the IAAPathComparator.
 * 
 * To be extended to implement alternative similarity variants
 * @author albertgo
 *
 */
public class HeadAtomComputer
{  
   public static final HeadAtomComputer INSTANCE = new HeadAtomComputer();
   
   protected HeadAtomComputer()
   {}
   
   protected double getHeadAtomSim(int atIdx1, IAAPathComputerInterface m1, 
                                   int atIdx2, IAAPathComputerInterface m2)
   {  int aType1 = m1.getAtomType(atIdx1);
      int aType2 = m2.getAtomType(atIdx2);
   
      if(aType1 == aType2) return 1D;
   
      return 0D;
   }
}
