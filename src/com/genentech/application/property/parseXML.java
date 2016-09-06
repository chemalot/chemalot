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
package com.genentech.application.property;
import java.io.*;
import java.util.*;
//import org.jdom.Attribute;
//import org.jdom.Comment;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.JDOMException;
import org.jdom.input.SAXBuilder;

/* parses an XML file and returns a LinkedHashMap with the
 * attribute name as key and attribute value as value*/
public class parseXML {

   /* Use jdom to read and parse an XML */
   private Document readDocument(String file) {
      try {
         SAXBuilder builder = new SAXBuilder();
         Document aDocument = builder.build(file);
         return aDocument;
      } catch(JDOMException e) {
         System.err.println(file + " is not well-formed.");
         System.err.println(e.getMessage());
         e.printStackTrace();
      } catch(NullPointerException e) {
         e.printStackTrace();
      } catch (IOException e) {
         System.err.println("Could not check " + file);
         System.err.println(" because " + e.getMessage());
      }
      return null;
   }

/* recursive procedure to add elements of the XML file to an hash */
   @SuppressWarnings("unchecked")
   private void addChildrenToHash(Element current, int depth, Map<String, String> elementsMap ) {

      //printSpaces(depth);
      //System.out.println(current.getName());
      if ( current.getName().compareTo("smarts") == 0) {
         elementsMap.put(current.getAttributeValue("name") , current.getTextTrim());
      }
      List<Element> children = current.getChildren();
      Iterator<Element> iterator = children.iterator();
      while (iterator.hasNext()) {
         Element child = iterator.next();
         addChildrenToHash(child, depth+1, elementsMap);
      }
   }

   @SuppressWarnings("unused")
   private void printSpaces(int n) {

      for (int i = 0; i < n; i++) {
         System.out.print(' ');
      }
   }

   public void printHash(Map<String, String> map) {
      Iterator<String> iterator = map.keySet().iterator();
      while(iterator.hasNext()){
         String name = iterator.next();
         String smarts = map.get(name);
         System.out.println(name + " " + smarts);
      }
   }


   /*function to parse a XML file and return a LinkedHashMap */
   /*map keys are attribute names and map values are attribute values */
   public LinkedHashMap<String, String> parse(Element root) {
      try {
         LinkedHashMap<String, String> elementsMap = new LinkedHashMap<String, String>();
         addChildrenToHash(root,0, elementsMap);
         //printHash(linkedHashMap);
         return (elementsMap);
      }
      catch (Exception e) {
         e.printStackTrace();
      }
      return null;
   }

   
   /*function to parse a XML file and return a LinkedHashMap */
   /*map keys are attribute names and map values are attribute values */
   public LinkedHashMap<String, String> parse(String fileName) {
      try {
         LinkedHashMap<String, String> elementsMap = new LinkedHashMap<String, String>();
         Document doc = readDocument(fileName);
         Element root = doc.getRootElement();
         addChildrenToHash(root,0, elementsMap);
         //printHash(linkedHashMap);
         return (elementsMap);
      }
      catch (Exception e) {
         e.printStackTrace();
      }
      return null;
   }

   public static void main(String args[]){
      if (args.length == 0) {
         System.err.println("Usage: java parseXML <fileName>");
         return;
      }

      parseXML myParseXML = new parseXML();
      try {
         Map<String, String> temp = myParseXML.parse(args[0]);
         myParseXML.printHash(temp);
      }
      catch (Exception e) {
         e.printStackTrace();
      }
   }
}
