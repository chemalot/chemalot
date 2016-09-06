#!/usr/bin/env perl
use warnings;
$script=$0;
$script =~ /^\// || ( $script="$ENV{PWD}/$script" );
our($installDir) = $script =~ /(.*)\/[^\/]+$/;

use Getopt::Long;
use File::Copy;
use File::Path;

#Example commands are in bin/R subfolder

$progName = "sdfRRandomForestCreator.pl\n";
$use = <<___help;
sdfRRandomForestCreator.pl
   -in <arg>   Supported file format: .sdf, .tab
   -out <arg>   Supported file format: .sdf, .tab
   -modelLocation <arg> If not specified, the model and the 
         companion files will be deposited in
         /gne/research/data/smdd/ApplicationData/ModelData. If you
         wish to store the output files in a sub directory specify
         the path to the desired sub directory. To save the files
         elsewhere, enter relative ("." or "..") or the absolute path.
   -modelName <arg>   The file extension .RData will be appended to
         the model file. Two companion text file will be generated as
         well. [modelFileName]_descriptor.txt contains the
         descriptorFields input as entered.
         [modelFileName]_stdout.log contains output from stdout.
   -identifierField <arg>   The field/column containing the record
         (compound) identifier. If is not given the creator will
         generated the identifiers.
   -responseField <arg>   The field/column containing the
         to-be-predicted values.
   -descriptorFields <arg>   Separate the name with |, e.g f1|f2.
   -predictionColName <arg>   Name of the output column with the
         prediction. Required if no conversion Script
   -nTree <arg>   Number of trees to grow. If not specified,
         default is 500
   -mTry <arg>   Number of variables randomly sampled as candidates
         at each split. If not specified the default for
         classification is sqrt(p) and for regression is p/3 where
         p is number of descriptors.
   -nodeSize size at which trees stop growing. Default 5.
   -logImportance   If specify, the importance of descriptors will
         be assessed and output to stdout.
   -modelOptions <arg>   Any RandomForest-specific arguments in the
         format as described in the R Documentation for "Classification
         and Regression with Random Forest separated by commas.
   -seed <arg>   integer seed for R.
   -correlationCutoff <arg>   Integer value used as the threshold
         to remove descriptors that correlate with other another
         descriptors above the given cutoff value. If no cutoff is
         specified, no descriptor will be removed. Note that
         descriptor with zero variance will always be removed.
   -conversionScript <arg>   File name plus the path to the script.
         This R script contains two functions. The function
         "convertInput" is used by this script to convert the input
         response before creating the model. The function 
         "descriptorManipulation" can be used to manipulate the descriptor
         matrix eg. by calculating quatratic terms. The function
         "computeOutput" is stored and to be used by the RModelPredictor
         to generate the additional output. If not specified, response
         is used as it is and no instruction for additional output is
         stored for the RModelScript to use.
   -printRScript   Output R script before execution.
   -writeInTab <fname>   To write the input table to a file for
         debugging
   -h   Output the help page
___help


#Following setting for near neigbor analysis
#   should be the same in all RModelCreator and RModelPredictor scripts
$IDENTIFIERTag = "_internalID_";
$SDFCfpLevel = 2;
$FINGERPRINTTag = "AFP".$SDFCfpLevel;


#For logging the command line executed in this script
$inputCommand = "sdfRRandomForestCreator.pl ".join( " ", @ARGV );

$modelSpecificOutputCommand = "";

$nTree = "";
$mTry = "";
$nodeSize = "5";
$logImportance = "";
$printRScript = "";
$predictionColName = "";
GetOptions( 'in=s' => \$input,
            'out=s' => \$output,
            'modelLocation=s' => \$modelLocation,
            'modelName=s' => \$modelName, 
            'identifierField=s' => \$identifier,
            'responseField=s' => \$response,
            'descriptorFields=s' => \$descriptor, 
            'predictionColName=s' => \$predictionColName,
            'nTree=i' => \$nTree, 
            'mTry=s' => \$mTry, 
            'nodeSize=s' => \$nodeSize, 
            'logImportance' => \$logImportance,
            'modelOptions=s' => \$modelOptions,
            'seed=s' => \$seed,
            'correlationCutoff=f' => \$correlationCutoff,
            'conversionScript=s' => \$conversionScript,
            'printRScript' => \$printRScript,
            'writeInTab=s' => \$writeInTab,
            'h' => \$help );
if( $help )
{  warn( "$use");
   exit();
}
warn( "$progName \n\n" );


$errMessage = "";
if( !$input )                                                                   
{  $errMessage .= "in is reqired\n"; }
if( !$modelLocation )
{  $errMessage .= "modelLocation is reqired\n"; }
if( !$modelName )
{  $errMessage .= "modelName is reqired\n"; }
if( !$response )
{  $errMessage .= "responseField is reqired\n"; }
if( !$descriptor )
{  $errMessage .= "descriptorFields is reqired\n"; }
if( !$predictionColName && !$conversionScript )
{  $errMessage .= "predictionColName or conversionScript is reqired\n"; }
if( $conversionScript && !-e $conversionScript )
{  $errMessage .= "$conversionScript does not exists\n";
}
$errMessage && die( "$errMessage \n$use" );


#Auxilary files for R computation
$sdfRCreatorFunctionsFile = "$installDir/R/model/sdfRCreatorFunctions.R";

umask 0007;
$modelRootDir = "/gne/research/data/smdd/ApplicationData/ModelData";
if( $modelLocation =~ m/^[^.\/]/ )
{  $modelDir = "$modelRootDir/$modelLocation";
} else
{  $modelDir = $modelLocation;
}
if( !-d $modelDir )
{  mkpath( $modelDir, {mode=>0770} ); 
}
$modelFile = "$modelDir/$modelName";

#Logging the commandline executed in this script
$commandLogFile = "$modelFile.commandLog.txt";
printToCommandLogfile( "$commandLogFile", "$inputCommand" );

$descriptorFile = "$modelFile.descriptor.txt";
$logRStdoutFile = "$modelFile.stdout.log";
$trainingFile   = "$modelDir/_trainingSet";

#Copy the conversion RScript to the model directory
#scriptFile: User's original file
$scriptCopy = "$modelFile.conversionFunctions.R";
if( $conversionScript )
{  system( "cp $conversionScript $scriptCopy" );
} else
{  $scriptCopy = "";
}


if( !$modelOptions )
{  $modelOptions = "";
}
if( $nTree )
{ $modelOptions .= ", ntree=$nTree";
}
if( $mTry )
{ $modelOptions .= ", mtry=$mTry";
}
if( $nodeSize )
{ $modelOptions .= ", nodesize=$nodeSize";
}
if( $logImportance )
{  $modelOptions .= ", importance=TRUE";
   $modelSpecificOutputCommand = "round( importance( model ), 2 )";
}

if( $modelOptions )
{  $modelOptions = ", $modelOptions";
}
#die( "modelOptions = $modelOptions \n" );

#R command for setting the seed in R
$setSeedCommand = "";
if( $seed )
{  $setSeedCommand = "set.seed( $seed )";
}

if( !$correlationCutoff )
{  $correlationCutoff = "-1"; }


#Name Pipe did not work; ToDo: Need more trials.
$tmpDir = "/tmp";
$tmpTrainingFile = "$tmpDir/_trainingSet.$$";
$dataFile = "$tmpDir/Rmc.$$";
if( $input =~ /\.tab$/i )
{  my $option = "-addCounter $IDENTIFIERTag";
   if( $identifier )
   {  $option = "-rename $identifier=$IDENTIFIERTag";
   }
   $trainingFile = "$trainingFile.tab";
   $tmpTrainingFile = "$tmpTrainingFile.tab";
   if( $input eq ".tab" )
   {  $com = "cat - | tabTagTool.pl $option > $trainingFile";
   } else
   {  -e $input || die "$input does not exists";
      $com = "cat $input | tabTagTool.pl $option > $trainingFile";
   }
   printToCommandLogfile( "$commandLogFile", "$com" );
   system( $com );
   system( "tabTagTool.pl $trainingFile -keep '$descriptor|$response' > $dataFile" );
   printToCommandLogfile( "$commandLogFile", "$com" );
   system( $com );
    
} elsif( $input =~ /\.sdf$/i )
{  $trainingFileTSMI = "$trainingFile.TSMI.tab";
   $trainingFile = "$trainingFile.sdf";
   $tmpTrainingFile = "$tmpTrainingFile.sdf";
   if( $input ne ".sdf" )
   {  -e $input || die "$input does not exists";
   }
   $addCounter_com = "";
   if( !$identifier )
   {  $identifier = "counter";
      $addCounter_com = " | sdfTagTool.csh -in .sdf -out .sdf -addCounter";
   }
   $setID_com = $addCounter_com
               ." | sdfTagTool.csh -in .sdf -out .sdf -rename '"
               .$identifier."=".$IDENTIFIERTag."'";
   $com = "sdfCFP.csh -in $input -out .sdf -level $SDFCfpLevel"
         .$setID_com
         ." | sdfNormalizer.csh -in .sdf -out .sdf"
         ." | tee $trainingFile"
         ." | sdf2Tab.csh -in .sdf > $trainingFileTSMI -tags 'CTISMILES|$IDENTIFIERTag'";
   printToCommandLogfile( "$commandLogFile", "$com" );
   system( $com );
   $com = "sdf2Tab.csh -in $trainingFile -tags '$descriptor|$response' >> $dataFile";
   printToCommandLogfile( "$commandLogFile", "$com" );
   system( $com );
} else
{  die( "Unknown file extension type" );
}


if( $writeInTab )
{  system( "cp $dataFile $writeInTab" );
}


# R will autoconvert hyphens to "."
my( $responseForR ) = $response;
$responseForR =~ s/-/./g;

$com = <<com;
   R --vanilla --slave >& "$logRStdoutFile"  <<'SCRIPT'
      writeLines( " " )
      writeLines( "Model name: $modelName" )
      writeLines( "Location: $modelDir" )
      
      #Contains functions required by this script for the R model creation
      source( "$sdfRCreatorFunctionsFile" )
      
      #Contains user defined the function for prepossing input response and 
      #the default function for generating the output after the prediction
      sourceConversionScript( "$scriptCopy" )
         
      dataFile <- "$dataFile"
      modelFile <- paste( "$modelFile", '.RDATA', sep='' )
      responseField <- "$responseForR"
      predictionColName <- "$predictionColName"
      logImportance <- $logImportance
      $setSeedCommand
      
      data <- read.table( dataFile, sep='\\t', header=TRUE, as.is=TRUE )

      # save RDATA so we can debug if necessary
      save.image(sub("RDATA", "savePoint.RDATA", modelFile))

      data <- removeNAandINF( data )
      
      y <- as.numeric(as.vector(data[, responseField]))
      #Convert the response, if specified
      if( exists( "convertInputResponse", mode="function" ) )
      {  y <- convertInputResponse( y )
         remove( convertInputResponse )
      }
      x <- data[, colnames( data ) != responseField]
      # transform descriptors if requested
      if( exists( "descriptorManipulation", mode="function" ) )
      {  x <- descriptorManipulation( x )
      }

      # save RDATA so we can debug if necessary
      save.image(sub("RDATA", "savePoint.RDATA", modelFile))
      
      correlationCutoff <- $correlationCutoff

      # save RDATA so we can debug if necessary
      save.image(sub("RDATA", "savePoint.RDATA", modelFile))

      x <- removeCorrelatedData( x, correlationCutoff )
      
      # preset variable p which can be used to define mTry on command line
      p <- ncol(x)
      
      #modelOptions is a perl string variable containing the model-specific 
      #option setting. The perl script has to validate the option specification.
      libList<-c( "randomForest" )
      library( "randomForest" )
      
      # save RDATA so we can debug if necessary
      save.image(sub("RDATA", "savePoint.RDATA", modelFile))
      
      # Remove to save memory for model creation
      remove(data)

      model <- randomForest( x, y $modelOptions )
      
      #Compute Rsquare (out-of-bag) and RMS error (out-of-bag)
      Rsquare <- cor( model\$predicted, y, use="complete.obs", method="pearson" )^2
      RMS <- sqrt( mean( ( model\$predicted - y )^2, na.rm = TRUE ) )
      cat( "\nR square (out-of-bag):", format( Rsquare, digits=4 ) )
      cat( "\nRMS error (out-of-bag):", format( RMS, digits=4 ), "\n" )
      
      model
      
      writeLines( sprintf( "\nSize of training set: %i", nrow(x) ) )
      writeLines( sprintf( "Number of descriptors: %i \n", ncol(x) ) )
      
      #The perl script should set, e.g. round( importance( model ), 2 ) for RandomForest
      $modelSpecificOutputCommand
      
      modelDescriptors <- colnames(x)
      save.image(sub("RDATA", "savePoint.RDATA", modelFile))
      remove( dataFile, responseField, x, y )
      save.image(modelFile)

      R.version.string
      quit( status=0 )
SCRIPT
com


if( $printRScript )
{  warn( $com );
}
$status = system( $com );

# print R stderr
open(INFILE, "$logRStdoutFile") || die "$!";
while($_=<INFILE>) 
{  warn $_; 
}
close(INFILE);

if( $status != 0 ) 
{  die "Error in R\n";
}

system( "echo '".$descriptor."' > ".$descriptorFile );


#Add the prediction of just created model
if( $trainingFile =~ /\.sdf$/i )
{  $com = <<___coms;
      sdfRModelPredictor.pl -in $trainingFile -out $tmpTrainingFile \\
            -modelLocation $modelLocation -modelName $modelName \\
            -FPFieldName $FINGERPRINTTag
___coms
   printToCommandLogfile( "$commandLogFile", "$com" );
   $status = system( $com );
   if( $status == 0 ) 
   {  $com = "mv $tmpTrainingFile $trainingFile";
      printToCommandLogfile( "$commandLogFile", "$com" );
      system( $com ); 
   }
}

# Output the training data if requested
if( $output )
{  if( $output =~ /^\./i )
   {  system( "cat $trainingFile" );
   } else
   {  system( "cp $trainingFile $output" );
   }
}



#Write critical command text to $commandLogFile
sub printToCommandLogfile
{  #Function variable declaration and assignment to the local variables.
   my( $commandLogFile, $LogText ) = @_;
   open( COMMANDLogFile, ">>$commandLogFile" ); 
   print( COMMANDLogFile "$LogText \n\n" );
   close( COMMANDLogFile );
}


