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


public enum MCSSCompareType
{  DEFAULT,
   QUERYRatio;
   
   public static MCSSCompareType toEnum( String str )
   {  String s = str.toUpperCase();
      for( MCSSCompareType enumvar : MCSSCompareType.values() )
         if( enumvar.name().toUpperCase().equals(s)) 
            return enumvar;
      
      throw new Error(String.format("Unknonw MCSSComapreType: %s", str));
   }
}
