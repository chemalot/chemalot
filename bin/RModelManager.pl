#!/usr/bin/env perl
use warnings;
use Getopt::Long;
use File::Copy;
use File::Path;

#Example commands are in ~smdi/dev/bin/R subfolder

$progName = "RModelManager.pl\n";
$use = <<___help;
$progName
   -modelLocation <arg>   If not specified, the model and the companion
         files is assumed to be in
         /gne/research/data/smdd/ApplicationData/ModelData/RModel. If
         the model output files are in a sub directory, specify the
         path to the given sub directory. If the model output files are
         elsewhere, enter relative ("." or "..") or the absolute path
   -modelName <arg>   The name of the model output files.
   -computeOutputScript <arg>   The R script containing the functions
         "computeOutput". The function "computeOutput" in this script
         will overwrite the existing function for generating the
         additional output. The specification is optional.
   -printRScript   Output R script before execution.
   -h   Output the help page
___help


$printRScript = "";
$debug        = "";
$sdfRCreatorFunctionsFile = "$ENV{AESTEL_DIR}/../bin/R/model/sdfRCreatorFunctions.R";

umask 0007;



GetOptions( 'modelLocation=s' => \$modelLocation,
            'modelName=s' => \$modelName, 
            'computeOutputScript=s' => \$computeOutputScript,
            'printRScript' => \$printRScript,
            'debug' => \$debug,
            'h' => \$help ) || die("$use\n");
if( $help )
{  warn( "$use");
   exit();
}
warn( "$progName \n\n" );


$errMessage = "";
if( !$modelLocation )
{  $errMessage .= "modelLocation is reqired\n"; }
if( !$modelName )
{  $errMessage .= "modelFilename is reqired\n"; }
$errMessage && die( "$errMessage \n$use" );


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

if( !$computeOutputScript )
{  $computeOutputScript = "";
}elsif( ! -e $computeOutputScript )
{  die "conversionScript does not exist: $computeOutputScript\n";
}

$errMessage = "";
if( !-e "$modelFile.RDATA" )
{  die( "$modelFile.RDATA required\n$use");
}
system("cp '$modelFile.RDATA' '$modelFile.RDATA.bak'") == 0 || die "$!";

$com = <<com;
   R --vanilla --slave --args  <<'SCRIPT'   
      source( "$sdfRCreatorFunctionsFile" )

      args <- commandArgs(TRUE)
      # argument index starts with 1 instead 0
      load( paste("$modelFile", '.RDATA', sep='') )
      
      if( length( "$computeOutputScript" ) > 0 
       && file.access( "$computeOutputScript", mode=4 ) == 0 )
      {  if( exists( "computeOutput", mode="function" ) )
            remove( computeOutput )
         if( exists( "getOutputNameText", mode="function" ) )
            remove( getOutputNameText )
         sourceConversionScript( "$computeOutputScript" )
      }
      remove(args)
      save.image( paste("$modelFile", '.RDATA', sep='') )

      quit(status=0)
SCRIPT
com


if( $printRScript || $debug )
{  warn( $com );
}
exec( $com );

