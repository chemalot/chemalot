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

import java.io.*;
import java.util.*;
import java.util.Map.Entry;
import java.util.regex.Pattern;

import com.aestel.Settings;
import com.aestel.chemistry.openEye.fp.tools.GeneralHashFunctionLibrary;
import com.aestel.utility.LogHelper;
import com.aestel.utility.Message;
import com.aestel.utility.Message.Level;
import com.aestel.utility.exception.AestelError;

/**
 * Use a tab separated file to map structure Code names to positions in a bitString.
 * If addNewCodes is false this is not mutable.
 *
 * The dictionary can contain bits of various types to create combined fingerprints.
 *
 * @author albertgo
 *
 */
public class TABDictionaryStructureCodeMapper implements LearningStrcutureCodeMapper
{  /**
    * Map from the type of a StructureCodeIterator to a Map which contains the
    * position of the bit keyed by the StructureCodeName.
    */
   private final Map<String,Map<String,Integer>> typeToNameMapMap
                                          = new HashMap<String, Map<String, Integer>>();
   private final String fileName;
   private final boolean addNewCodes;
   private final boolean hashNewCodes;
   private final List<String> indxToName = new ArrayList<String>(2000);
   private int maxIdx = -1;
   private int minIdx = 0;
   private int numHashIdexes = 0;


   /**
    * A {@link StructureCodeMapper} which derives the position from a dictionary
    * stored in an tab separated file in AESTEL_DIR/config/fp. Depending on the value
    * of addNewCodes codeNames which are not contained in the dictionary are
    * either ignored or appended to the list of known codeNames.
    *
    * The columns in the file are: type,idx,codeName
    * @param fileName filename of the tab separated dictionary inside the fp directory.
    * @param addNewCodes if false getPosition returns -1 for codes with unknown name.
    *
    */
   public TABDictionaryStructureCodeMapper(String fileName, boolean addNewCodes, boolean hashNewCodes)
   {  this.fileName = fileName;

      if(addNewCodes && hashNewCodes)
         throw new AestelError(new Message("addNewCodes and hashNewCodes may not be used together.",
                                           Level.ERROR, null));
      this.addNewCodes = addNewCodes;
      this.hashNewCodes = hashNewCodes;

      fileName = Settings.AESTEL_INSTALL_PATH + "/config/fp/" + fileName;
      File file = new File( fileName );
      if(! file.exists() )
      {  LogHelper.LOG.warning("Fingerprint file not found: " + file.toString());
         return;
      }

      minIdx = Integer.MAX_VALUE;
      maxIdx = Integer.MIN_VALUE;
      Pattern commentPat = Pattern.compile("//.*");
      Pattern tabPat = Pattern.compile("\t");
      try
      {  BufferedReader in = new BufferedReader(new FileReader(file));
         String line;
         while((line=in.readLine()) != null)
         {  line = commentPat.matcher(line).replaceFirst("");
            line = line.trim();
            if(line.length() == 0) continue;

            String[] vals = tabPat.split(line);
            if(vals.length != 3)
            {  LogHelper.LOG.warning(String.format("Invalid line in %s: %s",fileName,line));
               continue;
            }
            String type = vals[0];
            int idx = Integer.parseInt(vals[1]);
            String name = vals[2];

            Map<String,Integer> nameToIdxMap = typeToNameMapMap.get(type);
            if(nameToIdxMap == null)
            {  nameToIdxMap = new HashMap<String, Integer>();
               typeToNameMapMap.put(type, nameToIdxMap);
            }
            nameToIdxMap.put(name, idx);

            while( indxToName.size() <= idx )
               indxToName.add(null);
            if( indxToName.get(idx) != null && indxToName.get(idx).equals(name))
               System.err.printf("Index reused: %s %s\n", indxToName.get(idx), name);
            else
               indxToName.set(idx, name);

            if(idx > maxIdx) maxIdx = idx;
            if(idx < minIdx) minIdx = idx;
         }
         in.close();
      } catch (Exception e)
      {  LogHelper.LOG.severe( "FileName: " + file.toString() );
         e.printStackTrace();
         throw new Error(e);
      }
      if(minIdx == Integer.MAX_VALUE) minIdx = 0;
      if(maxIdx == Integer.MIN_VALUE) maxIdx = -1;

      if(hashNewCodes)
      {  if(maxIdx < 100)
            throw new Error("hashNewCodes should only be used with larger dictionaries\n");
         numHashIdexes = (maxIdx-minIdx)/3;
      }

//      for(Entry<String,Integer> e : typeToNameMapMap.get("lin74").entrySet())
//         System.err.println(e.getKey() + " " + e.getValue());
   }

   /**
    * returns -1 if addNewCodes is false and the codeName is unknown.
    */
   @Override
   public int getIndex(String type, String codeName)
   {  Map<String, Integer> nameToIdxMap = typeToNameMapMap.get(type);
      if(nameToIdxMap == null)
      {  // unknown fragment type
         if(addNewCodes)
         {  nameToIdxMap = new HashMap<String,Integer>();
            typeToNameMapMap.put(type, nameToIdxMap);
         }else if(hashNewCodes)
         {  int idx = GeneralHashFunctionLibrary.toInt(
               GeneralHashFunctionLibrary.BKDRHash(codeName), numHashIdexes)+maxIdx+1;
            return idx;
         }else
         {  return -1;
         }

      }

      // known fragment type
      Integer idx = nameToIdxMap.get(codeName);
      if( idx != null) return idx;

      // unknown fragmentCodename
      if(addNewCodes)
      {  nameToIdxMap.put(codeName, ++maxIdx);
         return maxIdx;
      } else if(hashNewCodes)
      {  idx = GeneralHashFunctionLibrary.toInt(
            GeneralHashFunctionLibrary.BKDRHash(codeName), numHashIdexes)+maxIdx+1;
         return idx;
      }else
      {  return -1;
      }
   }

   /**
    * Writes the dictionary into the dictionary file.
    * The original file is renamed to .bak.
    */
   @Override
   public void writeDictionary()
   {  String myfileName = Settings.AESTEL_INSTALL_PATH + "/config/fp/" + fileName;
      String destName = myfileName.replaceAll("\\.tab$", "") + ".bak";
      File file = new File( myfileName );
      File destFile = new File(destName);
      destFile.delete();
      file.renameTo(destFile);
      destName = myfileName;

      // sort by types
      String[] types = typeToNameMapMap.keySet().toArray(new String[typeToNameMapMap.size()]);
      Arrays.sort(types);


      Comparator<Entry<String, Integer>> compareByPos = new Comparator<Entry<String, Integer>>()
      {  @Override
      public int compare(Entry<String, Integer> a, Entry<String, Integer> b)
         { return a.getValue().compareTo(b.getValue());}
      };

      try
      {  PrintWriter out = new PrintWriter(new FileWriter(file));
         out.println("//type\tidx\tcodeName");

         for(String type : types)
         {  Map<String, Integer> nameToIdxMap = typeToNameMapMap.get(type);
            List<Entry<String, Integer>> nameIdxList
                  = new ArrayList<Entry<String, Integer>>(nameToIdxMap.size());
            nameIdxList.addAll(nameToIdxMap.entrySet());
            Collections.sort(nameIdxList, compareByPos);

            for( Entry<String, Integer> nameIdx : nameIdxList )
            {  String name = nameIdx.getKey();
               int pos = nameIdx.getValue();

               out.printf("%s\t%d\t%s\n", type,pos,name);
            }
         }
         out.close();

      } catch (Exception e)
      {  LogHelper.LOG.severe( "FileName: " + file.toString() );
         e.printStackTrace();
         throw new Error(e);
      }

   }

   @Override
   public void close()
   {  // nothing to do here
   }

   @Override
   public int getMaxIdx()
   {  return maxIdx + numHashIdexes;
   }

   @Override
   public int getMinIdx()
   {  return minIdx;
   }

   /**
    *
    * @return iterator of String[3] = {CodeType, CodeName, BitIndex}.
    */
   @Override
   public Iterator<String[]> getIterator()
   {  return new TABCodeIterator();
   }

   @Override
   public void reSortDictionary(int[] freq)
   {  class SCode { final String type; final String name; final int freq;
                    SCode(String type, String name, int freq)
                    {  this.type = type; this.name = name; this.freq = freq;
                    }
                    @Override
                  public String toString()
                    { return String.format("%s\t%s\t%d", type, name, freq); }
                  };

      SCode[] sCodes = new SCode[maxIdx+1];
      String[] types = typeToNameMapMap.keySet().toArray(new String[typeToNameMapMap.size()]);
      for (String type : types)
      {  Map<String, Integer> nameToIdxMap = typeToNameMapMap.get(type);
         for (Entry<String, Integer> nameIdx : nameToIdxMap.entrySet())
         {  int idx = nameIdx.getValue();
            SCode code = new SCode(type, nameIdx.getKey(), freq[idx]);
            sCodes[idx] = code;
         }
      }
      Comparator<SCode> cmp = new Comparator<SCode>()
      {  @Override
      public int compare(SCode o1, SCode o2)
         {  if( o1 == null) return 1;
            if( o2 == null) return -1;
            return o1.freq < o2.freq ? 1 : (o1.freq > o2.freq ? -1 : 0);
         }
      };
      Arrays.sort(sCodes, cmp);

      int idx = 0;
      for(int i=0; i<sCodes.length; i++)
      {  SCode code = sCodes[i];
         System.err.println(code);
         if(code == null) continue;
         Map<String, Integer> nameToIdxMap = typeToNameMapMap.get(code.type);
         nameToIdxMap.put(code.name, idx++);
      }
   }

   /**
    * iterator over all codes in this dictionary as Array of 3 strings (type, name, bitIndex).
    * @author albertgo
    *
    */
   class TABCodeIterator implements Iterator<String[]>
   {  private Iterator<String> typeIterator;
      private Iterator<Entry<String, Integer>> codeIterator;
      private String[] code = new String[3];
      public int hashIdx = 0;

      public TABCodeIterator()
      {  typeIterator = typeToNameMapMap.keySet().iterator();

         assert typeIterator.hasNext();
         code[TYPEIdx] = typeIterator.next();
         codeIterator = typeToNameMapMap.get(code[TYPEIdx]).entrySet().iterator();
      }
      @Override
      public boolean hasNext()
      {  if(code[INDEXIdx] != null) return true;
         if(codeIterator.hasNext()) {
            Entry<String, Integer> tmp = codeIterator.next();
            code[NAMEIdx] = tmp.getKey();
            code[INDEXIdx] = tmp.getValue().toString();
            return true;
         }

         if(typeIterator.hasNext())
         {  code[TYPEIdx] = typeIterator.next();
            codeIterator = typeToNameMapMap.get(code[TYPEIdx]).entrySet().iterator();
            return hasNext();
         }

         if(hashIdx < numHashIdexes)
         {  code[INDEXIdx] = Integer.toString(maxIdx+(++hashIdx));
            code[NAMEIdx] = "HASH" + code[INDEXIdx];
            code[TYPEIdx] = "HASHED";
            return true;
         }

         return false;
      }


      @Override
      public String[] next()
      {  if(!hasNext())
            throw new NoSuchElementException();

         String[] tmp = code;
         code = tmp.clone();
         code[INDEXIdx] = null;

         return tmp;
      }

      @Override
      public void remove()
      {  throw new UnsupportedOperationException();
      }
   }

   @Override
   public String getName(int idx)
   {  assert idx < indxToName.size();

      return indxToName.get(idx);
   }
}
