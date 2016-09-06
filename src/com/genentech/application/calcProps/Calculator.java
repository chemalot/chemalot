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
package com.genentech.application.calcProps;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.jdom.Element;

/**
 * This class stores information about a property calculator needed for generation
 * of a command line text.
 * 
 * Native calculators are generated using definitions in properties files. Some
 * of the calculators contains aggregation ID. Calculators with the same aggregation
 * ID can be consolidated into an aggregation calculator because the associated 
 * properties can be computed by the same program instance.
 *
 * @author JW Feng & Man-Ling Lee / Last updated on Oct 31, 2015
 * Copyright 2012-2015 Genentech
*/
public class Calculator
{
   private String name=null;
   private boolean isPublic=false;
   private boolean requiresIonization=false;
   private String progAggregateID;
   private String progName=null;
   private String progOps=null;
   private Set<String> requiredCalculators = new HashSet<String>();
   private Set<String> keepRequiredCalculators = new HashSet<String>();
   private List<String> outputFields = new ArrayList<String>();
   private List<String> verboseFields = new ArrayList<String>();
   private String helpText = null;
   private List<String> originalCalculatorList = new ArrayList<String>();
   

   public Calculator(Element propertyNode, boolean useExp)
   {
      String temp ="";
      Set<String> emptySet = new HashSet<String>( 0 );
      
      this.name = propertyNode.getAttributeValue("name");
      this.isPublic= "Y".equalsIgnoreCase(propertyNode.getAttributeValue("isPublic"));
      this.requiresIonization= "Y".equalsIgnoreCase(propertyNode.getAttributeValue("requiresIonization"));
      this.progAggregateID= propertyNode.getAttributeValue("progAggregateID") == null ? 
               "" :  propertyNode.getAttributeValue("progAggregateID").trim().replaceAll("\\s", "").toLowerCase().intern();
      
      String expProgName=null;
      String expProgOptions=null;
      expProgName = propertyNode.getChildTextTrim("expProgName");
      expProgOptions = propertyNode.getChildTextTrim("expProgOptions"); 
      if (expProgName == null && expProgOptions != null) { // must specify expProgName first
         System.err.println("Error in specifying options for " + this.name + " in the xml files.");
         System.err.println("Cannot specify an expProgOption without specifying a expProgName.");
         System.exit(1);
      }

      String tempProgName=null;
      String tempProgOps=null;
      tempProgName = propertyNode.getChildTextTrim("progName");
      tempProgOps = propertyNode.getChildTextTrim("progOptions"); 
      if (tempProgName == null && tempProgOps != null) { // must specify progName first
         System.err.println("Error in specifying options for " + this.name + " in the xml files.");
         System.err.println("Cannot specify an progOption without specifying a progName.");
         System.exit(1);
      }
      
      //by default use progName and progOps
      this.progName= tempProgName == null ? "" : tempProgName.trim();
      this.progOps= tempProgOps == null ? "" : tempProgOps.trim();
      
      if (useExp && expProgName != null) 
      { // expProgName if available and useExp is specified
         this.progName = expProgName == null ? "" : expProgName.trim();
         this.progOps= expProgOptions == null ? "" : expProgOptions.trim();
      }
      
      temp = propertyNode.getChildTextTrim("requiredCalculators");
      if (temp != null && temp.length() > 0) {
         Collections.addAll( this.requiredCalculators, temp.trim().split("\\|") );  

         if (this.requiresIonization == true) {
           System.err.println("Calculators that require ionization cannot have requiredCalculator.");
           System.err.println("Please fix entry for: " + this.name);
           System.exit(1);
         }
      }
      
      temp = propertyNode.getChildTextTrim("keepRequiredCalculators");
      if (temp != null && temp.length() > 0)
         Collections.addAll( this.keepRequiredCalculators, temp.trim().split("\\|") );
      else
         this.keepRequiredCalculators = emptySet;
         
      //if requiredCalculators is not defined, use keepRequiredCalculators
      if( this.requiredCalculators == null || this.requiredCalculators.size() == 0 )
         this.requiredCalculators.addAll( this.keepRequiredCalculators );
      
      temp=propertyNode.getChildTextTrim("outputFields");
      if (temp != null && temp.length() > 0)
         Collections.addAll( this.outputFields, temp.trim().split("\\|") );
         
      temp=propertyNode.getChildTextTrim("verboseFields");
      if (temp != null && temp.length() > 0)
         Collections.addAll( this.verboseFields, temp.trim().split("\\|") );
         
      temp = propertyNode.getChildTextTrim("helpText");
      this.helpText = temp == null ? "" : temp.trim();
   }

   
   /**
    * @return the name
    */
   public String getName()
   {
      return name;
   }
   void setName( String propertyName )
   {
      name = propertyName;
   }
   /**
    * @return the progName
    */
   public String getProgName()
   {
      return progName;
   }
   /**
    * @return the progNameOps
    */
   public String getProgOps()
   {
      return progOps;
   }
   
   /**
    * @param ops
    * @return 
    */
   public void setProgOps(String ops){
      this.progOps= ops;
   }
   
   /**
    * @return the reqiredCalculators
    */
   public Set<String> getRequiredCalculators()
   {  return requiredCalculators; }
   
   public void setRequiredCalculators( Set<String> calculatorNames )
   {  this.requiredCalculators = calculatorNames; }
   
   void addRequiredCalculators( Set<String> calculatorNames )
   {  requiredCalculators.addAll( calculatorNames ); }
  
   /**
    * @return the keepRequiredCalculators
    */
   public Set<String> getKeepRequiredCalculators()
   {  return keepRequiredCalculators; }
   
   void addKeepRequiredCalculators( Set<String> calculatorNames )
   {  keepRequiredCalculators.addAll( calculatorNames ); }
   

   /**
    * @return the isPublic
    */
   public boolean isPublic()
   {  return isPublic; }

   
   /**
    * @return the outputFields
    */
   public List<String> getOutputFields()
   {  return outputFields; }
   
   void addOutputFields( List<String> fields )
   {  outputFields.addAll( fields ); }

   /**
    * @return the verboseOutputFields
    */
   public List<String> getVerboseFields()
   {  return verboseFields; }
   
   void addVerboseFields( List<String> fields )
   {  verboseFields.addAll( fields ); }

   
   /**
    * @return the progAggregateID
    */
   public String getProgAggregateID()
   {  return progAggregateID; }


   /**
    * @return the helpText
    */
   public String getHelpText()
   {  return helpText; }
   
   void setHelpText( String text )
   {  this.helpText = text; }
   
   
   List<String> getOriginalCalculators()
   {  return this.originalCalculatorList; }
   
   void addOriginalCalculator( String calculatorName )
   {  this.originalCalculatorList.add( calculatorName ); }
   

   /**
    * @return the requiresIonization
    */
   public boolean requiresIonization()
   {  return requiresIonization; }

 


}
