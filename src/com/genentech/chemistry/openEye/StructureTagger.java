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
package com.genentech.chemistry.openEye;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;

import openeye.oechem.OEMatchBaseIter;
import openeye.oechem.OEMolBase;
import openeye.oechem.OESubSearch;



/**
 * Use it to tag a structure based on SMARTS or molfile input.
 * 
 * After creating the StructureTagger, structure patterns for tagging have to be 
 * added first. Call tagStructure() to execute the tagging, then call the get 
 * methods to retrieve the tagging information. Before the next execution of 
 * the tagging, call the clear() method to remove the existing tagging information.
 *
 * @author Man-Ling Lee / February 11, 2011
 * Copyright 2011-2012 Genentech
 */
public class StructureTagger
{  private List<StructurePattern> patternList;
   /**
    * List of structure patterns that are found in the compound structure
    */
   private List<StructurePattern> matchedPatterns;
  
   
   public StructureTagger()
   {  this.patternList     = new ArrayList<StructurePattern>();
      this.matchedPatterns = new ArrayList<StructurePattern>();
   }
   
   
   public void clear()
   {  patternList.clear();
      matchedPatterns.clear();
   }
   
   
   public void addPattern( String patternString, OESubSearch subSearch, 
            String tagName, String setName )
   throws NullPointerException
   {  StructurePattern pattern = 
            new StructurePattern( patternString, subSearch, tagName, setName );
      patternList.add( pattern );
   }
   
   
   public boolean tagStructure( OEMolBase mol )
   {  matchedPatterns.clear();
      boolean matched = false;
      for( int i=0; i<patternList.size(); i++ )
      {  OEMatchBaseIter match = patternList.get(i).subSearch.Match( mol, true );
         if( match.hasNext() )
         {  matchedPatterns.add( patternList.get(i) );
            matched = true;
         }
         match.delete();
      }
      return matched;
   }
   
   public LinkedHashMap<String,Integer> countOccurrence( OEMolBase mol )
   {  LinkedHashMap<String,Integer> map = new LinkedHashMap<String,Integer>();
      for( int i=0; i<patternList.size(); i++ )
      {  OEMatchBaseIter match = patternList.get(i).subSearch.Match( mol, true );
         int count = 0;
         while( match.hasNext() )
         {  match.next().delete();
            ++count;
         }
         match.delete();
         String tagName = patternList.get(i).tagName;
         Integer newCount = map.get( tagName );
         if( newCount == null )
            newCount = count;
         else
            newCount = newCount + count;
         map.put( patternList.get(i).tagName, newCount );
      }
      return map;
   }
   
   
   public LinkedHashMap<String,Boolean> checkOccurrence( OEMolBase mol )
   {  LinkedHashMap<String,Boolean> map = new LinkedHashMap<String,Boolean>();
      for( int i=0; i<patternList.size(); i++ )
      {  boolean doesMatch = patternList.get(i).subSearch.SingleMatch(mol);
         String tagName = patternList.get(i).tagName;
         
         Boolean oldVal = map.get( tagName );
         if( oldVal == null )
            map.put(tagName, Boolean.valueOf(doesMatch));
         else if( doesMatch && oldVal.booleanValue() == false )
            map.put(tagName, Boolean.valueOf( Boolean.TRUE ));
      }
      return map;
   }
   
   
   public int getTagCounts()
   {  return matchedPatterns.size();
   }
   
   
   public String getAllTags()
   {  if( matchedPatterns.size() == 0 )
         return "";
      StringBuffer sb = new StringBuffer();
      for( int i=0; i<matchedPatterns.size(); i++ )
         sb.append( matchedPatterns.get(i).tagName ).append( ';' );
      return sb.substring( 0, sb.length()-1 );
   }
   
   
   public String getAllTagPattern()
   {  if( matchedPatterns.size() == 0 )
         return "";
      StringBuffer sb = new StringBuffer();
      for( int i=0; i<matchedPatterns.size(); i++ )
         sb.append( matchedPatterns.get(i).patternString ).append( ' ' );
      return sb.substring( 0, sb.length()-1 );
   }
   
   
   public String getFirstTag()
   {  if( matchedPatterns.size() == 0 )
         return "";
      return matchedPatterns.get(0).tagName;
   }
   
   
   public String getFirstTagPattern()
   {  if( matchedPatterns.size() == 0 )
         return "";
         return matchedPatterns.get(0).patternString;
   }

   

   static class StructurePattern
   {  final String patternString;
      final String tagName;
      final String setName;
      final OESubSearch subSearch;
   
      StructurePattern( String patternString, OESubSearch subSearch, 
               String tagName, String setName )
      throws NullPointerException
      {  if( subSearch == null || tagName == null || tagName.length() == 0 )
            throw new NullPointerException( "Missing pattern or tagName input." );
         if( patternString == null )
            this.patternString = "";
         else
            this.patternString = patternString;
         this.subSearch = subSearch;
         this.tagName = tagName;
         if( setName == null )
            this.setName = "";
         else
            this.setName = setName;
      }
   }
}
