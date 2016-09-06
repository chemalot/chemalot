# Do not change the name of the three functions!!!!

#Will be call by the R model creator to generate the input for the model creation. 
convertInputResponse <- function( responseVector )
{  return( log10(responseVector) + 3 )
}


# Will be call by the R model predictor after the prediction to generate the output
# fields
#
# predictedDataMatrix is assumed to have one column with the predicted values 
# with the rownames(predictedDataMatrix) containing the row identifier.
#
# This method should output a data table containing computed columns
computeOutput <- function( predictedDataMatrixi, data )
{  #Compute probability of being less than the threshold, if lower.tail = TRUE,
   # i.e. integration from minus infinity to threshold
   exp_stddev <- 0.5
   exp_threshold <- 4
   probability <- pnorm( predictedDataMatrix, mean = exp_threshold, sd = exp_stddev, 
                         lower.tail = TRUE )
   
   
   #Compute output in addition to probability from the predicted values
   solubility_mM <- 10^( predictedDataMatrix - 6 )
   cLogS <- predictedDataMatrix - 9
   
   
   #Compute categorical output derived from the continuous output
   solubility_cat <- vector()
   solubility_cat[ solubility_mM <= 10 ] <- 'low'
   solubility_cat[ solubility_mM <= 100 ] <- 'medium'
   solubility_cat[ solubility_mM > 100 ] <- 'high'
   
   
   #Return the computed values and the row names copied from the input matrix
   df<-data.frame( sol_prob=probability, 
                   solubility_uM=solubility_uM,
                   solubility_cat,
                   cLogS=cLogS )
   
   
   #Don't change the two lines below
   rownames( df ) <- rownames( predictedDataMatrix )
   return( df )
}


# Will be call by the R model predictor after the prediction to retrieve the names
# of the corresponding computed columns
#
# This method should output a String with tab-delimited column names
getOutputNameText <- function()
{  return( "sol_prob\tsolubility_uM\tsolubility_cat\tcLogS" )
}

