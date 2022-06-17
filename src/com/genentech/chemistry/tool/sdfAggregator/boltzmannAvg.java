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
package com.genentech.chemistry.tool.sdfAggregator;

import java.util.ArrayList;
import java.util.Collection;

/**
 * This class implements the mean function.
 * @author Ben Sellers 2017-05-23
 *
 */

import com.aestel.utility.DataFormat;

import openeye.oechem.OEGraphMol;
import openeye.oechem.oechem;

public class boltzmannAvg extends AggFunction
{  
   // Members for energy calculations needed for Bolztmann weighting
   protected String energyTag;
   protected Collection<Double> energyContainer;
   protected double minEnergy = Double.MAX_VALUE;
   protected final double RT = 0.5921856D;  // at 298K
 
   
   public static final String DESCRIPTION
         = "[outName[:format] = ] boltzmannAvg(propertyTag,energyTag): returns the Boltzmann-weighted average of 'propertyTag' given the energy in 'energytag'.\n";


   public boltzmannAvg(String outTag, String resultFormat, String funcName, String funcArg)
   {  
      super(outTag, resultFormat == null ? "si2" : resultFormat, funcName, funcArg);    
      this.energyContainer = new ArrayList<Double>();       
   }

   @Override
   protected String getTagName(String funcArg)
   {  
      // Expecting two tag names, one for the property and one for the energy
      funcArg = funcArg.trim();
      if (!funcArg.contains(","))
         return "";
      
      String propertyTag = funcArg.substring(0, funcArg.indexOf(",")).trim();
      this.energyTag     = funcArg.substring(funcArg.indexOf(",")+1).trim();
            
      return propertyTag;
   }

   private String getResult(String fmtStr)
   {  
      if( valueContainer.size() == 0 )
         return "";
      
      if( energyContainer.size() == 0 )
         return "";
      
      // Must be energy for every value/property
      if( energyContainer.size() != valueContainer.size() )
         return "";      
              
      double numeratorSum = 0D;
      double deltaE = Double.MAX_VALUE; 
      double energy = Double.MAX_VALUE;
      double bFactor = Double.MAX_VALUE;
      double partitionFunc = 0;
      int valueIndex = 0; 
   
      
      for(String v : valueContainer)
      {  
         // Get energies
         try
         {  
            energy = (double) energyContainer.toArray()[valueIndex];
            deltaE = energy - this.minEnergy;
            bFactor = Math.exp(-1 * deltaE / RT);
            partitionFunc+=bFactor; 
         } 
         catch (NumberFormatException e)
         {  
            System.err.printf("Warning: %s cannot be converted to float.\n", e);
         }
         
         try
         {  
            numeratorSum += (Float.parseFloat(v) * bFactor);
         } 
         catch (NumberFormatException e)
         {  
            System.err.printf("Warning: %s for %s cannot be converted to float.\n",
                              v, getOutTagName());
         }
         valueIndex++;
      }
      
      double expectationValue = numeratorSum / partitionFunc;

      return DataFormat.formatNumber(expectationValue, fmtStr);
   }
  

   @Override
   public String getResult(int indxInGrp)
   {  return getResult(resultFormat);
   }
   
   @Override
   public void process (OEGraphMol mol)
   {  
      String propertyTagVal = oechem.OEGetSDData(mol, this.tagName);
      if( propertyTagVal != null && propertyTagVal.length() > 0)
         valueContainer.add(propertyTagVal);
      
      // Also get the energy
      String energyValStr = oechem.OEGetSDData(mol, this.energyTag);
      if( energyValStr != null && energyValStr.length() > 0)
      {
         // Parse it double
         double energyVal = Double.MAX_VALUE; 
         try
         {  
             energyVal = Double.parseDouble(energyValStr);
             energyContainer.add(energyVal);
             checkMinEnergy(energyVal);
         } 
         catch (NumberFormatException e)
         {  
            System.err.printf("Warning: %s cannot be converted to float.\n", energyVal);
         }                 
      }
   }

   @Override
   public void init()
   {  
      valueContainer.clear();
      energyContainer.clear();
      minEnergy = Double.MAX_VALUE;
   }
   
   private void checkMinEnergy(Double energyVal)
   {      
       if (energyVal < this.minEnergy)
          this.minEnergy = energyVal;
           
   }
}
