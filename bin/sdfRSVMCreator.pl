#!/usr/bin/env perl
use warnings;
$script=$0;
$script =~ /^\// || ( $script="$ENV{PWD}/$script" );
our($installDir) = $script =~ /(.*)\/[^\/]+$/;

use Getopt::Long;
use File::Copy;
use File::Path;

#Example commands are in R subfolder

$progName = "sdfRSVMCreator.pl";
$use = <<___help;
$progName
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
         prediction. Required if no conversionScript.
   -modelOptions <arg>   Any SVM-specific arguments in the format
         as described in the R Documentation for "Classification
         and Regression with SVM separated by commas.
   -gamma <arg>   Parameter needed for all kernels except linear
         (default: 1/(data dimension)) Parameter neded for all kernel
         types other than Linear. Typically, this parameter must be
         tuned in order to get the best model. The default is one
         over the number of descriptor properties, nx. A single
         value, a comma-separated list of values, or any R expression
         generating a vector of values may be specified. If multiple
         values are given, a model will be built for each value. If
         CrossValidate is set to True, cross-validation will be
         performed for each value, and the value giving the best
         cross-validation score will be used to build the final model
         from all the data. In specifying values, you may use nx to
         indicate the number of descriptor properties.
   -cost <arg>   Cost of constraints violation (default: 1)it is the
         C-constant of the regularization term in the Lagrange
         formulation. Cost associated with training set errors.
         Typically, this parameter must be tuned in order to get
         the best model. A single value, a comma-separated list of
         values, or any R expression generating a vector of values
         may be specified. If multiple values are given, a model
         will be built for each value. If CrossValidate is set to
         True, cross-validation will be performed for each value,
         and the value giving the best cross-validation score will
         be used to build the final model from all the data.
   -crossValidationFolds <arg>   The number of 'folds' for 
         cross-validation. (default: 10) E.g., setting this parameter
         to 10 implies 10-fold cross-validation, in which 9/10 of
         the data are used to train a model and the remaining 1/10
         to test the model. The process is repeated 10 times so that
         each sample is left out once.
   -type <arg>   SVM can be used as a classification machine,
         as a regression machine, or for novelty detection. Depending
         of whether y is a factor or not, the default setting for
         type is C-classification or eps-regression, respectively,
         but may be overwritten by setting an explicit value. Valid
         options are: C-classification, nu-classification, 
         one-classification (for novelty detection), eps-regression,
         and nu-regression
   -epsilon <arg>   Epsilon in the insensitive-loss function
         (default: 0.1)
   -nu <arg>   Parameter needed for nu-classification, nu-regression,
         and one-classification.
   -kernel <arg>   Kernel function to be used for mapping input space
         to feature space. Valid: linear, polynomial, radial, and 
         sigmoid
   -degree <arg>   Parameter needed for kernel of type polynomial
         (default: 3)
   -coef0 <arg>   Parameter needed for kernels of type polynomial
         and sigmoid (default: 0)
   -seed <arg>   Integer seed for libsvm (used for cross-validation
         and probability prediction models)
   -correlationCutoff <arg>   Number value used as the threshold
         to remove descriptors that correlate with other another
         descriptors above the given cutoff value. If no cutoff is
         specified, no descriptor will be removed. Note that descriptor
         with zero variance will always be removed.
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
$inputCommand = "sdfRSVMCreator.pl ".join( " ", @ARGV );


$modelSpecificOutputCommand = "";
$tuneOptions = "";
$predictionColName = "";

$printRScript = "";
GetOptions( 'in=s' => \$input, 
            'out=s' => \$output,
            'modelLocation=s' => \$modelLocation,
            'modelName=s' => \$modelName,
            'identifierField=s' => \$identifier,
            'responseField=s' => \$response,
            'descriptorFields=s' => \$descriptor, 
            'predictionColName=s' => \$predictionColName,
            'modelOptions=s' => \$modelOptions,
            'gamma=s' => \$gamma,
            'cost=s' => \$cost,
            'crossValidationFolds=s' => \$crossValidationFolds,
            'epsilon=s' => \$epsilon,
            'nu=s' => \$nu,
            'kernel=s' => \$kernel,
            'degree=s' => \$degree,
            'coef0=s' => \$coef0,
            'seed=i' => \$seed,
            'type=s' => \$type,
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
{  $errMessage .= "predictionColName is reqired\n"; }
if( $conversionScript && !-e $conversionScript )
{  $errMessage .= "$conversionScript does not exists\n"; }
if( !$cost )                                                                   
{  $errMessage .= "cost is reqired\n"; }

#ToDo Instead, we could allow default types, i.e. C-classification or eps-regression,
#     depending of whether y is a factor or not.
if( !$type )
{  $errMessage .= "type is reqired\n"; }

if( $kernel && $kernel !~ m/linear|polynomial|radial|sigmoid/i )
{  $errMessage .= "Input kernel expression is wrong.\n"
                 ."Valid: linear, polynomial, radial, and sigmoid\n";
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


if( $gamma )
{  if( $gamma !~ m/\(/ )
   {  $gamma = "c($gamma)";
   }
   $tuneOptions .= ", gamma=$gamma";
}
if( $cost )
{  if( $cost !~ m/\(/ )
   {  $cost = "c($cost)";
   }
   $tuneOptions .= ", cost=$cost";
}
if( !$crossValidationFolds )
{  $crossValidationFolds = 5
}

if( !$modelOptions )
{  $modelOptions = "";
}
if( $epsilon )
{  $modelOptions .= ", epsilon=$epsilon";
}
if( $nu )
{  $modelOptions .= ", nu=$nu";
}
if( $kernel )
{  $kernel = lc( $kernel );
   $modelOptions .= ", kernel=\"$kernel\"";
}
if( $degree )
{  $modelOptions .= ", degree=$degree";
}
if( $coef0 )
{  $modelOptions .= ", coef0=$coef0";
}
if( $seed )
{  $modelOptions .= ", seed=$seed";
}

if( $modelOptions && $modelOptions !~ m/^,/ )
{  $modelOptions = ", $modelOptions";
}
#die( "modelOptions = $modelOptions \n" );

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
      
      correlationCutoff <- $correlationCutoff

      # save RDATA so we can debug if necessary
      save.image(sub("RDATA", "savePoint.RDATA", modelFile))

      x <- removeCorrelatedData( x, correlationCutoff )
      
      nx <- length(x)
      gamma <- $gamma
      cost <- $cost
      
      #modelOptions is a perl string variable containing the model-specific 
      #option setting. The perl script has to validate the option specification.
      libList<-c("e1071")
      library("e1071")
      
      # save RDATA so we can debug if necessary
      save.image(sub("RDATA", "savePoint.RDATA", modelFile))
      
      # Remove to save memory for model creation
      remove(data)

      if( length( gamma ) > 1 || length( cost ) > 1 )
      {  #Cross-validation to identify the best model
         tuneControl <- tune.control( sampling = "cross", 
               cross = $crossValidationFolds, best.model=TRUE )
         tuneResult <- tune( svm, x, y, 
               ranges = list( gamma = gamma, cost = cost ),
               tunecontrol = tuneControl $modelOptions )
         print( summary( tuneResult ) )
         
         #Get the model trained on all training data using the best parameter set
         model <- tuneResult\$best.model
         
         save.image(sub("RDATA", "savePoint.RDATA", modelFile))
         remove( tuneResult )
      } else 
      {  #Just create a model
         model <- svm( x, y 
               $modelOptions $tuneOptions, cross=$crossValidationFolds )
         save.image(sub("RDATA", "savePoint.RDATA", modelFile))
      }
      model
      
      writeLines( sprintf( "Size of training set: %i", nrow(x) ) )
      writeLines( sprintf( "Number of descriptors: %i \n", ncol(x) ) )
      
      #The perl script should set, e.g. "round( importance( model ), 2 )" for svm
      $modelSpecificOutputCommand
      
      modelDescriptors <- colnames(x)
      remove( dataFile, responseField, x, y )
      save.image( modelFile )

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
      sdfRModelPredictor.pl -in $trainingFile -out  $tmpTrainingFile \\
            -modelLocation $modelLocation --modelName $modelName \\
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


