# Remove alphanumeric columns and zero-variance columns as wells as columns with 
# data that correlated with data in another column from inData
removeCorrelatedData <- function( inData, correlationCutoff )
{
   dataColNames <- names( inData )
   
   #Keep the name of alphanumeric columns 
   #Create a vector containing the variance of each column
   #tmp <- diag( var( inData, na.rm=FALSE ) )
   tmp <- apply( inData, 2, var )
   txtColNames <- dataColNames[is.na( tmp )]
   
   #Remove the alpha numeric columns
   inData <- inData[!dataColNames %in% txtColNames]
   dataColNames <- names( inData )
   
   #Keep the names of columns with zero variance
   tmp <- as.data.frame( rbind( tmp[!is.na( tmp )] ) )
   tmpCols <- apply( tmp == 0, 2, any )
   zeroVarColNames <- names( tmpCols[tmpCols] )
   
   #Remove the zero variance columns
   inData <- inData[!dataColNames %in% zeroVarColNames]
   dataColNames <- names( inData )
   
   #Keep the names of columns that are correlated to an existing column
   if( correlationCutoff > -1 )
   {  tmp <- cor( data.frame( lapply( inData, as.numeric ), check.names=FALSE ), 
                              use="pairwise.complete.obs" )
      diag( tmp ) <- 0
      tmp[lower.tri( tmp )] <- 0
      tmpCols <- apply( abs( tmp ) >= correlationCutoff, 2, any )
      colinearColNames <- names( tmpCols[tmpCols] )
   
      #Remove the correclated columns
      inData <- inData[!dataColNames %in% colinearColNames]
      dataColNames <- names( inData )
   
      # Find names of removed columns and their corresponding most correlated column
      if( length( colinearColNames ) > 0 )
      {  pairsOfCorrelatedCols <- apply( 
               abs(tmp), 2, function(x) names( x[rank( x,TRUE,"random" ) == max( rank( x, TRUE, "first" ) )] ) )[tmpCols]
   
         #Output summary of removed columns
         writeLines( "\nColumns removed due to strong correlation with data in other columns:")
         print( data.frame(removed=names(pairsOfCorrelatedCols),correlated=pairsOfCorrelatedCols, row.names=1) )
         writeLines( "\n" )
      }
   }
   return( inData )
}



sourceConversionScript <- function( scriptFile )
{  if( length( scriptFile ) > 0 && file.access( scriptFile, mode=4 ) == 0 )
   {  source( scriptFile )
      if( ( exists( "computeOutput", mode="function" ) 
            && !exists( "getOutputNameText", mode="function" ) )
       || ( !exists( "computeOutput", mode="function" ) 
            && exists( "getOutputNameText", mode="function" ) ) )
      {  write( paste( "\nError:",
                       "Both functions computeOutput() and getOutputName()",
                       "have to be defined in", scriptFile, "or none at all.\n" ), 
                stderr() );
         quit( status=1 )
      }
   } 
}



anyNA <- function( inData ) 
{  return( any( is.na( inData ) ) ) }



#Issue warnings for following cases
#  1) are there any is.na or is.infinite values in the response column
#  2) are there any is.na or is.inifinite values in the descriptor columns
#
#  Remove rows by applying to the whole data set before split in to x and y
removeNAandINF <- function( inData )
{  dataColNames <- colnames( inData )
   dataRowNames <- rownames( inData )            
   
   tmpData <- as.matrix( inData )
   tmpData[ is.infinite( tmpData ) ] <- NA
   
   naColNames <- dataColNames[ apply( tmpData, 2, anyNA ) ]
   naRowNames <- dataRowNames[ apply( tmpData, 1, anyNA ) ]
   
   if( length( naColNames ) > 0 || length( naRowNames ) > 0 )
   {  write( "\n########################################", stdout() )
      write( "Result of NA and INFINITE values removal", stdout() )
      if( length( naColNames ) > 0 )
      {  write( "Following columns contain NA and INFINITE value:", stdout() )
         write( naColNames, stdout() );
      }
      if( length( naRowNames ) > 0 )
      {  write( paste( "Following rows contain NA and INFINITE value",
                       "and have been removed:" ),
                stdout() );
         write( naRowNames, stdout() );
      }
      write("########################################\n", stdout() )
   }
   inData <- na.omit( inData )
   return( inData )
}


