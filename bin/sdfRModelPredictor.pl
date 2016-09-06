#!/usr/bin/env perl
use warnings;
$script=$0;
$script =~ /^\// || ( $script="$ENV{PWD}/$script" );
our($installDir) = $script =~ /(.*)\/[^\/]+$/;

use Getopt::Long;
use File::Copy;
use File::Path;

#Example commands are in ~smdi/dev/bin/R subfolder

$progName = "sdfRModelPredictor.pl\n";
$use = <<___help;
$progName
   -in <arg>   Supported file format: .sdf
   -out <arg>   Supported file format: .sdf
   -modelLocation <arg>   If not specified, the model and the companion
         files is assumed to be in
         /gne/research/data/smdd/ApplicationData/ModelData/RModel. If
         the model output files are in a sub directory, specify the
         path to the given sub directory. If the model output files are
         elsewhere, enter relative ("." or "..") or the absolute path
   -modelName <arg>   The name of the model output files.
   -NNSimAbove <arg>   Tanimoto coefficient needed for retrieving
         compounds in the training set above the specified similarity.
         Default values is 0.75
   -NNFieldPrefix <arg>   Prefix for names of the fields containing the
         output from the near neigbor analysis. If no prefix is
         specified, no near neigbor analysis is performed.
   -FPFieldName <arg>   Name of the field containing the hex-encoded fingerprints
         for nearneigbor analysis. If not specified, the AFP2 (Genentech atom 
         fingerprint) is used. This field will only be used if FPFieldName is 
         specified.
   -computeOutputScript <arg>   The R script containing the functions
         "computeOutput". The function "computeOutput" in this script
         will overwrite the existing function for generating the
         additional output. The specification is optional.
   -printRScript   Output R script before execution.
   -h   Output the help page
   -debug Print additional info and do not delete temporary files.
___help


#Following setting for near neigbor analysis
#   should be the same in all RModelCreator and RModelPredictor scripts
$IDENTIFIERTag = "_internalID_";
$SDFCfpLevel = 2;
$FINGERPRINTTag = "AFP".$SDFCfpLevel;


#The default environment for user SDMI is prd.
#If you want access accdev2, issue "source ~smdi/bin/setDev.csh" before 
# running this script.

$printRScript = "";
$debug        = "";
$fpFieldName  = "";

$tmpDir = "/tmp";
$user   = $ENV{'USER'};

#Logging the commandline executed in this script
$commandLogFile = "$tmpDir/RcmdLog.$user.$$.txt";

$inputCommand = "sdfRModelPredictor.pl ".join( " ", @ARGV );
printToCommandLogfile( "$commandLogFile", "$inputCommand" );


GetOptions( 'in=s' => \$input, 
            'out=s' => \$output, 
            'modelLocation=s' => \$modelLocation,
            'modelName=s' => \$modelName, 
            'NNSimAbove=f' => \$nnSimilarAbove,
            'NNFieldPrefix=s' => \$nnFieldPrefix,
            'FPFieldName=s' => \$fpFieldName,
            'computeOutputScript=s' => \$computeOutputScript,
            'printRScript' => \$printRScript,
            'debug' => \$debug,
            'h' => \$help );
if( $help )
{  warn( "$use");
   exit();
}
warn( "$progName \n\n" );


$errMessage = "";
if( !$input )                                                                   
{  $errMessage .= "in is reqired\n"; }
if( !$output )                                                                   
{  $errMessage .= "out is reqired\n";
} elsif( $output !~ /.*\.sdf(\.gz)?/ )
{   $errMessage .= "output type must be .sdf or .sdf.gz\n";
}
if( !$modelLocation )
{  $errMessage .= "modelLocation is reqired\n"; }
if( !$modelName )
{  $errMessage .= "modelFilename is reqired\n"; }
$errMessage && die( "$errMessage \n$use" );


if( !$nnSimilarAbove )
{  $nnSimilarAbove = 0.75; }


$modelRootDir = "/gne/research/data/smdd/ApplicationData/ModelData";
if( $modelLocation =~ m/^[^.\/]/ )
{  $modelDir = "$modelRootDir/$modelLocation";
} else
{  $modelDir = $modelLocation;
}
if( !-d $modelDir )
{  die( "$modelDir does not exists. \n$use" );
}
$modelFile = "$modelDir/$modelName";
$descriptorFile = "${modelFile}.descriptor.txt";
$trainingFile   = "$modelDir/_trainingSet.sdf";

if( !$computeOutputScript )
{  $computeOutputScript = "";
}

$errMessage = "";
if( !-e "$modelFile.RDATA" || !-e $descriptorFile )
{  die( "$modelFile.RDATA and $descriptorFile are required for prediction\n$use");
}
open( FILE, $descriptorFile ) or die $!;
$descriptor = <FILE>;
chomp $descriptor;
close( FILE );
$nDescriptor = 0+(@_=split(/\|/, $descriptor));

$inputFile= "$tmpDir/Rin.$user.$$.sdf";
$dataFile = "$tmpDir/Rmc.$user.$$";

# make named pipe to start sdf2tab in background so we can run R in parrallel
-e $dataFile && unlink($dataFile);
system("mknod $dataFile p") == 0 || die "Could not create named pipe: $!";

#This is needed to combine the input data with the R prediction output
$MERGEFieldName = "counter";

if( $input =~ /\.sdf$/i )
{  if( $input ne ".sdf" && ! -e $input )
   {  die "$input does not exists";
   }
   # the perl line  below is just a security measure against invalid sdf records and
   # can be removed if those are not expected (201309)
   $com = <<___coms;
      sdfTagTool.csh -in $input -out .sdf -addCounter \\
        | tee $inputFile \\
        | sdf2Tab.csh -in .sdf -tags '$MERGEFieldName|$descriptor' \\
        | perl -ne '\@_=split(/\\t/); if(\$#_== $nDescriptor) {print} else{ warn "Invalid line: \$_"}' \\
        >> $dataFile &
___coms
   printToCommandLogfile( "$commandLogFile", "$com" );
   system( $com ) == 0 || die "$!";
} else
{  die( "Unknown file extension type" );
}

# if input is to stdin we need to store it into a tmp file
#   so that it can be used in merging with the predictions(tmpOutput)
# tmpOutput contains the tab separated fields: id\tPredictions
#  and can be writen to a named pipe for merging with the input sdf
$tmpOutput   = "$tmpDir/Rout.$user.$$";
$tmpSdfRMLog = "$tmpDir/sdfRM.$user.$$.log";
$com = <<com;
   R --vanilla --slave --args "$dataFile" "$tmpOutput" "$modelFile" >& "$tmpSdfRMLog" <<'SCRIPT'   
      doPrediction <- function(data)
      {  # transform descriptors if requested
         if( exists( "descriptorManipulation", mode="function" ) )
         {  data <- descriptorManipulation( data )
         }

         x <- data.matrix(data[, modelDescriptors])
         if( ncol(x) != length(modelDescriptors) )
         {  write(paste("expected\\n",modelDescriptors,"\\n\\nfound:",colnames(x)),
                  stderr())
            q(status=1)
         }

         p<-predict(model,x)
         p<-as.matrix(p)
         if(length(p) != nrow(data))
         {  write(paste("Expected to get as many predictions as inputs.",
                        "inputs=",nrow(data),"predictions=",nrow(p)),
                  stderr())
            q(status=1)
         }
         rownames(p) <- rownames(data)
         return(p)
      }


      ##################################### main script 

      args <- commandArgs(TRUE)
      # argument index starts with 1 instead 0
      modelFile <- paste(args[3], '.RDATA', sep='')
      
      #cd ~manle/ApplicationDev/DMPK_Model_Infrastructure/test_sdfRModel_cmdProg
      #dataFile<-"Rmc.14701"; output<-"t.tab"; modelFile<-"test_RF_train.RDATA"
      load(modelFile)
      for( l in libList ) library( l, character.only = TRUE)
      # assumptions: model            contains model object, 
      #              modelDescriptors contains descriptor names
      
      if( length( "$computeOutputScript" ) > 0 
       && file.access( "$computeOutputScript", mode=4 ) == 0 )
      {  if( exists( "computeOutput", mode="function" ) )
            remove( computeOutput )
         if( exists( "getOutputNameText", mode="function" ) )
            remove( getOutputNameText )
         sourceConversionScript( "$computeOutputScript" )
      }
      if( !exists( "computeOutput", mode="function" ) )
      {  computeOutput <- function( predictedDataMatrix, data)
         {  return( predictedDataMatrix )
         }
      }
      # make backwards compatible for computeOutput implementations that take only one argument
      if( length(formals(computeOutput)) != 2 )
      {  orgComputeOutput <- computeOutput
         computeOutput <- function( predictedDataMatrix, data )
         {  return( orgComputeOutput( predictedDataMatrix ) )
         }
      }

      if( exists( "getOutputNameText", mode="function" ) )
      {  predictionColName <- getOutputNameText()
      }
      
      dataFile <- args[1]
      output <- args[2]
      inCon <- fifo(dataFile, open='rt', blocking = TRUE)
      outCon <- file(output)
      open(outCon,"w")
      write(paste("$MERGEFieldName\\t",predictionColName,sep=""),outCon)
      data <- read.table(inCon, sep='\\t', header=1, as.is=TRUE, row.names=1, nrows=100 )
      orgColNames <- colnames(data)
      while ( nrow(data) > 0 )
      {  colnames(data) <- orgColNames
         complete <- complete.cases(data);
         if( ! all(complete) )
         {  mis <- data[!complete,]
            write(paste("Found rows with missing values in: ",
                        colnames(mis)[apply(is.na(mis),2,any)]), stderr());
            write(paste("Row Names with missing values: ", rownames(mis)), stderr());
            data <- data[complete,];
            if( nrow(data) == 0 ) next;
         }

         #save.image('t.R')
         p <- doPrediction(data)
         p <- computeOutput( p, data )
         write.table(p,sep='\\t',quote=FALSE,file=outCon,col.names=FALSE)
         #write(colnames(p), stderr())

         #Try-error need because if input record number is multiple of 100
         #the read.table call following the last record will cause an error because
         #there is no more record.
         data <- try(read.table(inCon, sep='\\t', header=0, as.is=TRUE, row.names=1, nrows=100 ))
         if(inherits(data, 'try-error')) data <- data.frame()
      } 
      close(outCon)
      close(inCon)

      quit(status=0)
SCRIPT
com


if( $printRScript || $debug )
{  warn( $com );
}
printToCommandLogfile( "$commandLogFile", "$com" );
$status = system( $com );

# print R stderr
open( INFILE, "$tmpSdfRMLog" ) || die "$!";
while( $_=<INFILE> ) 
{  warn $_; 
}
close( INFILE );

if( $status != 0 )
{  &exitWithInfo( "Error in R\n" ); 
}


#Create the output directory if needed (stopped here)
my $index = rindex( $output,  '/');
if( $index > 0 )
{  my $outputDir = substr( $output, 0, $index );
   if( !-d $outputDir )
   {  mkpath( $outputDir, {mode=>0770} ); 
   }
}

#Assume that if training set file is an SDFile, it contains AFP2 fingerprints
$nnSimLabel = substr( $nnSimilarAbove,2,2 );
if( length( $nnSimLabel ) < 2 )
{  $nnSimLabel .= "0";
}
$nnCommand = "";
$removeOption = "-remove $MERGEFieldName";
$renameOption = "";
if( -e "$trainingFile" && $nnFieldPrefix )
{  if( "$fpFieldName" )
   {  if(  "$fpFieldName" ne $FINGERPRINTTag )
      {  die( "Fingerprint in the training set differs from the input fingerprint"
             ."Required fingerprint is $FINGERPRINTTag.\n")
      }
      $fpFieldCom = "";
      $removeOption = "-remove '$MERGEFieldName|NNIdx'";
      
   } else #Need to compute the FP for the input compounds
   {  $fpFieldCom = " | sdfCFP.csh -in .sdf  -out .sdf -level $SDFCfpLevel";
      $removeOption = "-remove '$MERGEFieldName|NNIdx|$FINGERPRINTTag.'";
   }
   $nnCommand = $fpFieldCom
               ." | sdfFPNNFinder.csh -in .sdf -out .sdf"
               ." -fpTag $FINGERPRINTTag -idTag $IDENTIFIERTag"
               ." -maxNeighbors 1 -minSimilarity 0 -ref $trainingFile"
               ." -countSimilarAbove $nnSimilarAbove";
   $renameOption = "-rename 'NNSim=${nnFieldPrefix}_nearTc"
                  ."|NNId=${nnFieldPrefix}_nearGn"
                  ."|NNCount_${nnSimilarAbove}=$nnFieldPrefix"
                  ."_Tc${nnSimLabel}NN'";
}

$com = <<com;
sdfTabMerger.csh -sdf $inputFile -tab $tmpOutput -out .sdf \\
   -addEmptyValues -mergeTag $MERGEFieldName -mergeCol $MERGEFieldName \\
   $nnCommand \\
   | sdfTagTool.csh -in .sdf -out $output $removeOption \\
         $renameOption
com

printToCommandLogfile( "$commandLogFile", "$com" );
system( $com ) == 0 || &exitWithInfo("$!");


if( ! $debug ) 
{   unlink( $inputFile, $dataFile, $tmpOutput, $commandLogFile, $tmpSdfRMLog  );
} else
{  &exitWithInfo("");
}



sub exitWithInfo
{  #Function variable declaration and assignment to the local variables.
   my( $info ) = @_;
   die  "\n$info\n\n"
       ."inputFile=$inputFile\n"
       ."dataFile=$dataFile\n"
       ."tmpOutput=$tmpOutput\n"
       ."HOST=$ENV{HOST}\n"
       ."command logfile=$commandLogFile\n";
}



#Write critical command text to $commandLogFile
sub printToCommandLogfile
{  #Function variable declaration and assignment to the local variables.
   my( $commandLogFile, $LogText ) = @_;
   open( COMMANDLogFile, ">>$commandLogFile" ); 
   print( COMMANDLogFile "$LogText \n\n" );
   close( COMMANDLogFile );
}

