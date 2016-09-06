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

import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import com.genentech.chemistry.openEye.AAPathComparatorFact.AAPathCompareType;

public class AAPathComparator3Test extends AbstractAAPathComparatorTest
{  /** matrix with expected similarities of {@link AbstractAAPathComparatorTest#smis}
    *  Lower Triangle will be copied in setUp.
    */
   protected double[][] myExpectedSims =
         {  {     1D,     0D,     0D,     0D,     0D,     0D,     0D,     0D,     0D,     0D },
            {     0D,     1D,     0D, .1111D, .0625D,     0D,     0D,     0D, .0067D,     0D },
            {     0D,     0D,     1D,     0D, .0625D,     0D,     0D,     0D, .0077D,     0D },
            {     0D,     0D,     0D,     1D, .3125D,     0D,     0D,     0D, .0328D,     0D },
            {     0D,     0D,     0D,     0D,     1D,     0D,     0D,     0D, .1117D,     0D },
            {     0D,     0D,     0D,     0D,     0D,     1D, .1282D, .2051D, .0540D, .2423D },
            {     0D,     0D,     0D,     0D,     0D,     0D,     1D, .3205D, .3990D, .2423D },
            {     0D,     0D,     0D,     0D,     0D,     0D,     0D,     1D, .1106D, .0751D },
            {     0D,     0D,     0D,     0D,     0D,     0D,     0D,     0D,     1D, .2423D },
            {     0D,     0D,     0D,     0D,     0D,     0D,     0D,     0D,     0D,     1D },
         };

   @BeforeClass
   public void setUp()
   {  super.setUp(myExpectedSims);
   }

   @Override
   @Test
   public void testComparator()
   {  super.testComparator();
   }

   @Override
   @AfterClass
   public void close()
   {  super.close();
   }


   @Override
   protected AAPathComparatorFact getComparatorFact()
   {  AAPathComparatorFact cFact = new AAPathComparatorFact(
                  AAPathCompareType.DEFAULT, 3);
      return cFact;
   }
}
