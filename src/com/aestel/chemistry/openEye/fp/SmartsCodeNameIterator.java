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

import java.io.File;
import java.util.*;

import org.jdom.Document;
import org.jdom.Element;
import org.jdom.input.SAXBuilder;

import openeye.oechem.OEMolBase;
import openeye.oechem.OESubSearch;

import com.aestel.Settings;
import com.aestel.utility.LogHelper;
import com.aestel.utility.Message;
import com.aestel.utility.Message.Level;
import com.aestel.utility.exception.AestelError;

/**
 * Iterates over all names of smarts matching an input molecule (@link init).
 *
 * The list of smarts is defined by the xml file.
 * @author albertgo
 *
 */
public class SmartsCodeNameIterator implements StructureCodeNameIterator
{  private final OESubSearch[] matcher;
   private final String[]      name;
   final String        type;
   private int nextMatchIdx;
   private String matchedCodeName;
   private OEMolBase mol;

   /**
    *
    * @param smartsXML defining smarts used to identify structure features.
    */
   @SuppressWarnings("unchecked")
   public static StructureCodeNameIterator createFromXML(String smartsXML)
   {  String fileName = Settings.AESTEL_INSTALL_PATH + "/config/fp/" + smartsXML;
      File file = new File( fileName );
      if(! file.exists() )
      {  throw new Error("Smarts file not found: " + file.toString());
      }

      String type = "";
      List<String> nameList   = new ArrayList<String>();
      List<String> smartsList = new ArrayList<String>();
      try
      {  SAXBuilder builder = new SAXBuilder(false);  //non validating

         Document xmlDoc = builder.build(file);
         Element root = xmlDoc.getRootElement();
         type = root.getAttributeValue("type");
         assert type != null && type.length()>0 : "Type is not defined";
         for(Element smart : (List<Element>)root.getChildren("smartsDef"))
         {  String name   = smart.getAttributeValue("name");
            String smarts = smart.getAttributeValue("smarts");
            nameList.add(name);
            smartsList.add(smarts);
         }
      } catch (Exception e)
      {  LogHelper.LOG.severe( "FileName: " + file.toString() );
         e.printStackTrace();
         throw new Error(e);
      }

      return new SmartsCodeNameIterator( type,
             nameList.toArray(new String[nameList.size()]),
             smartsList.toArray(new String[smartsList.size()]));
   }



   private SmartsCodeNameIterator(String type, String[] names, String[] smarts)
   {  this.type = type;
      name = names;

      int i=0;
      matcher = new OESubSearch[name.length];
      for(String smart : smarts)
      {  matcher[i] = new OESubSearch(smart);
         if(! matcher[i].IsValid() )
            throw new AestelError(
                     new Message("Invalid smarts: " + smart,Level.ERROR, null));
         i++;
      }

      assert name.length == i;
   }

   @Override
   public void init(OEMolBase mol)
   {  this.mol = mol;
      nextMatchIdx = 0;
      matchedCodeName = null;
   }

   @Override
   public boolean hasNext()
   {  if( matchedCodeName != null ) return true;

      while(nextMatchIdx < matcher.length)
      {  if(matcher[nextMatchIdx].SingleMatch(mol))
         {  matchedCodeName =  name[nextMatchIdx++];
            return true;
         }
         nextMatchIdx++;
      }

      return false;
   }

   @Override
   public String next()
   {  if(! hasNext() )
         throw new NoSuchElementException();

      String tmp = matchedCodeName;
      matchedCodeName = null;
      return tmp;
   }

   @Override
   public void close()
   {  for(OESubSearch sub : matcher)
         sub.delete();
   }

   /**
    * @throws UnsupportedOperationException;
    */
   @Override
   public void remove()
   {  throw new UnsupportedOperationException();
   }

   @Override
   public String getType()
   {  return type;
   }
}
