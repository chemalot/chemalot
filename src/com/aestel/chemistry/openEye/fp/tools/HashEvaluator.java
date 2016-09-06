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
package com.aestel.chemistry.openEye.fp.tools;

import java.io.*;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

public class HashEvaluator
{  static final long MAPSize = 16348;
   private static long NOSIGN = 0x7fffffffffffffffL;

   public static void main(String...args) throws IOException
   { System.out.println("hash\tsec\tNBuckets\tmaxCollision\tnCollisionBuckets\tavgCollisions\tworstBuket");
     for(int i=0; i<=12; i++)
         tryhash(i);
   }
   
   private static final void tryhash(int type) throws IOException
   {  ArrayList<Set<String>> hashKeys = new ArrayList<Set<String>>((int) MAPSize);
      BufferedReader in = new BufferedReader(new FileReader("c:\\home\\albertgo\\Aestel\\config\\fp\\linear74frag.txt"));
      long start = System.currentTimeMillis();
      
      for(int i=0; i<MAPSize; i++)
         hashKeys.add(new HashSet<String>(20));
      
      String s;
      while((s=in.readLine()) != null)
      {  long hash = doHash(type, s);
         int idx = (int)((hash & NOSIGN) % MAPSize);
         Set<String> keys = hashKeys.get(idx);
         keys.add(s);
      }
      in.close();
      double time = (System.currentTimeMillis()-start)/1000D;
      
      int maxCollisions = 0;
      int sumCollisions = 0;
      int maxCollisionIdx = -1;
      int nBuckets = 0;
      int nCollisionBuckets = 0;
      for(int i=0; i<hashKeys.size(); i++)
      {  Set<String> keys = hashKeys.get(i);
         int nCollisions = keys.size();
         if(nCollisions > 0) 
         {  nBuckets++;
            nCollisions--;
         }
         if(nCollisions > 0) nCollisionBuckets++;
         
         if(nCollisions > maxCollisions)
         {  maxCollisions = nCollisions;
            maxCollisionIdx = i;
         }
         sumCollisions += nCollisions;
      }
      
      System.out.printf("%s\t%.3f\t%d\t%d\t%d\t%f\t%d\t",
            hashName(type),time,
            nBuckets, maxCollisions, nCollisionBuckets, 
            sumCollisions/(float)MAPSize, maxCollisionIdx);
      System.out.println(hashKeys.get(maxCollisionIdx));
   }
   
   public static long doHash(int type, String s)
   {  switch(type)
      {  case 0:   return GeneralHashFunctionLibrary.APHash(s);
         case 1:   return GeneralHashFunctionLibrary.BKDRHash(s);
         case 2:   return GeneralHashFunctionLibrary.BPHash(s);
         case 3:   return GeneralHashFunctionLibrary.DEKHash(s);
         case 4:   return GeneralHashFunctionLibrary.DJBHash(s);
         case 5:   return GeneralHashFunctionLibrary.ELFHash(s);
         case 6:   return GeneralHashFunctionLibrary.FNVHash(s);
         case 7:   return GeneralHashFunctionLibrary.JSHash(s);
         case 8:   return GeneralHashFunctionLibrary.PJWHash(s);
         case 9:   return GeneralHashFunctionLibrary.RSHash(s);
         case 10:  return GeneralHashFunctionLibrary.SDBMHash(s);
         case 11:  return s.hashCode();
         case 12:  return new Hash(s).hash_code;
      }
      throw new Error();
   }

   public static String hashName(int type)
   {  switch(type)
      {  case 0:   return "APHash";
         case 1:   return "BKDRHash";
         case 2:   return "BPHash";
         case 3:   return "DEKHash";
         case 4:   return "DJBHash";
         case 5:   return "ELFHash";
         case 6:   return "FNVHash";
         case 7:   return "JSHash";
         case 8:   return "PJWHash";
         case 9:   return "RSHash";
         case 10:  return "SDBMHash";
         case 11:  return "java";
         case 12:  return "hashClass";
      }
      throw new Error();
   }
}
