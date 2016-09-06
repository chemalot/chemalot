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

public class IAAPathComparatorCharTest extends AbstractAAPathComparatorTest
{  /** matrix with expected similarities of {@link AbstractAAPathComparatorTest#smis}
    *  Lower Triangle will be copied in setUp.
    */
   protected double[][] myExpectedSims =
         {  {     1D,     0D,     0D,     0D,     0D,     0D,     0D,     0D,     0D,     0D },
            {     0D,     1D,     0D, .0588D, .0323D,     0D,     0D,     0D, .0033D,     0D },
            {     0D,     0D,     1D,     0D, .0323D,     0D,     0D,     0D, .0039D,     0D },
            {     0D,     0D,     0D,     1D, .1126D,     0D,     0D,     0D, .0088D,     0D },
            {     0D,     0D,     0D,     0D,     1D,     0D,     0D,     0D, .0336D,     0D },
            {     0D,     0D,     0D,     0D,     0D,     1D, .0373D, .0645D, .0148D, .0869D },
            {     0D,     0D,     0D,     0D,     0D,     0D,     1D, .1101D, .1767D, .0869D },
            {     0D,     0D,     0D,     0D,     0D,     0D,     0D,     1D, .0328D, .0211D },
            {     0D,     0D,     0D,     0D,     0D,     0D,     0D,     0D,     1D, .0869D },
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
                  AAPathCompareType.DEFAULT, 8);
      return cFact;
   }
}
