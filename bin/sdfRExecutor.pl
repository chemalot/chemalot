#!/usr/bin/env perl
use warnings;
use strict;

my($script)=$0;
$script =~ /^\// || ( $script="$ENV{PWD}/$script" );
our($installDir) = $script =~ /(.*)\/[^\/]+$/;

use Getopt::Long;
use File::Copy;
use File::Path;

#Example commands are in ~smdi/dev/bin/R subfolder

my($progName) = "sdfRExecutor.pl\n";
my($use) = <<___help;
$progName
   -in <arg>      Supported file format: .sdf
   -out <arg>     Supported file format: .sdf
   -RScript <arg> RScript to be executed. path is checked relative to current 
                  and relative to $installDir/R/RExecutor 
                  For example look at $installDir/R/RExecutor/R2.sdfR
   -dataFields    Pipe separated list of datafields to pass to R.
   -chunkSize n   If given the input is passed in chunks to R to save memory
   -removeQualifier If given "[ ><~=]" will be removed from the input fields
   -printRScript   Output R script before execution.
   -params str    A "," separated list of additional parameter to pass to docompute
                  method. e.g. "method='spearman',use='all.obs'"
   -h   Output the help page
   -debug Print additional info and do not delete temporary files.
___help


my($printRScript)    = "";
my($debug)           = "";
my($removeQualifier) = "";
my($chunkSize)       = "";

my($tmpDir)= "/tmp";
my($user)  = $ENV{'USER'};

#Logging the commandline executed in this script
my($commandLogFile)= "$tmpDir/RcmdLog.$user.$$.txt";

my($inputCommand)= "sdfRExecutor.pl ".join( " ", @ARGV );
printToCommandLogfile( "$commandLogFile", "$inputCommand" );

my($rScript) = "";
my($input) = "";
my($output) = "";
my($dataFields) = "";
my($computeOutputScript) = "";
my($params) = "";
my($help) = "";
if($#ARGV >= 0 && $ARGV[0] =~ '.sdfR$') 
{  $rScript=shift;
}

GetOptions( 'in=s' => \$input, 
            'out=s' => \$output, 
            'RScript=s' => \$rScript,
            'dataFields=s' => \$dataFields, 
            'removeQualifier' => \$removeQualifier,
            'computeOutputScript=s' => \$computeOutputScript,
            'chunkSize=i' => \$chunkSize,
            'params=s' => \$params,
            'printRScript' => \$printRScript,
            'debug' => \$debug,
            'h' => \$help );

my($errMessage) = "";
if( !$rScript )
{  $errMessage .= "rScript is reqired\n"; 
}
my($dummy) = "$installDir/R/RExecutor/$rScript";
if( ! -e $rScript && ! -e $dummy)
{  $errMessage .= "rScript ($rScript) not found\n"; 
} else
{  -e $rScript || ($rScript = $dummy);
   $rScript = &readFile($rScript);
}
if( $help )
{  if( ! $errMessage )
   {  $rScript =~ s/^#\!.*$(\s+)//m;  # remove hashbang line
      $rScript =~ s/^[^#].*//sm; # keep only initial block of comments
      $rScript =~ s/^# ?//gm;
      die( "\n$rScript\n");
   }
      
   &exitWithHelp($use);
}
warn( "$progName\n\n" );


if( !$input )                                                                   
{  $errMessage .= "in is reqired\n"; 
}
if( !$output )                                                                   
{  $errMessage .= "out is reqired\n";
}
if( !$dataFields )                                                                   
{  $errMessage .= "dataFields is reqired\n";
}
my($nDataFields) = split(/\|/,$dataFields);
$errMessage && &exitWithHelp( $use, $errMessage );

if($params) { $params = ",$params "; }

if( !$computeOutputScript )
{  $computeOutputScript = "";
}

my($inputFile)= "$tmpDir/Rin.$user.$$.sdf";
my($dataFile) = "$tmpDir/Rmc.$user.$$";
$chunkSize && ($chunkSize = ",nrows=$chunkSize");

# make named pipe to start sdf2tab in background so we can run R in parrallel
-e $dataFile && unlink($dataFile);
system("mknod $dataFile p") == 0 || die "Could not create named pipe: $!";

#This is needed to combine the input data with the R output
my($MERGEFieldName) = "counter";

if( $input !~ /^\./ && ! -e $input )
{  die "$input does not exists";
}
if( $removeQualifier ) { $removeQualifier = "perl -pe 's/ *[<>=~] *//g'"; }
my($com) = <<___coms;
   sdfTagTool.csh -in $input -out .sdf -addCounter \\
   | tee $inputFile \\
   | sdf2Tab.csh -in .sdf -tags '$MERGEFieldName|$dataFields' \\
   $removeQualifier \\
   | perl -ne '\@_=split(/\\t/); if(\$#_!= $nDataFields) {print} else{ warn "Invalid line: \$_"}' \\
   >> $dataFile &
___coms
printToCommandLogfile( "$commandLogFile", "$com" );
system( $com ) == 0 || die "$!";

# if input is to stdin we need to store it into a tmp file
#   so that it can be used in merging with the predictions(tmpOutput)
# tmpOutput contains the tab separated fields: id\tPredictions
#  and can be writen to a named pipe for merging with the input sdf
my($tmpOutput)   = "$tmpDir/Rout.$user.$$";
my($tmpSdfRMLog) = "$tmpDir/sdfRM.$user.$$.log";
# TODO once hmeasure is installed in global lib add --vanilla to R options
$com = <<com;
   R --slave --args "$dataFile" "$tmpOutput" >"$tmpSdfRMLog" 2>&1 <<'SCRIPT'
      
      ##################################### main script 
      args <- commandArgs(TRUE)
      # argument index starts with 1 instead 0

      dataFile <- args[1]
      output <- args[2]
      inCon <- fifo(dataFile, open='rt', blocking = TRUE)
      outCon <- file(output)
      open(outCon,"w")

      ###### inserted below user supplied doCompute function
      $rScript
      ###### inserted above user supplied doCompute function
   
      # R write.table does nto output header for rowNames
      cat("$MERGEFieldName\\t",file=outCon)

      data <- read.table(inCon, sep='\\t', header=1, as.is=TRUE, row.names=1 $chunkSize )
      #save.image('t.RDATA')
      orgColNames <- colnames(data)
      printColNames=TRUE
      while ( nrow(data) > 0 )
      {  colnames(data) <- orgColNames
         
         #save.image('t.R')
         p <- doCompute(data $params)
         write.table(p,sep='\\t',quote=FALSE,file=outCon,col.names=printColNames)
         printColNames = FALSE
         
         #Try-error need because if input record number is multiple of chunkSize
         #the read.table call following the last record will cause an error because
         #there is no more record.
         data <- try(read.table(inCon, sep='\\t', header=0, as.is=TRUE, row.names=1 $chunkSize ))
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
my($status) = system( $com );

# print R stderr
open( INFILE, "$tmpSdfRMLog" ) || die "$!";
while( $_=<INFILE> ) 
{  warn $_; 
}
close( INFILE );

if( $status != 0 )
{  &exitWithInfo( "Error in R\n" ); 
}


$com = <<com;
sdfTabMerger.csh -sdf $inputFile -tab $tmpOutput -out .sdf \\
   -addEmptyValues -mergeTag "$MERGEFieldName" -mergeCol "$MERGEFieldName" \\
   | sdfTagTool.csh -in .sdf -out $output -remove "$MERGEFieldName"
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

sub exitWithHelp
{  my($help, $msg) = @_;
   $msg || ($msg = "");
   $msg && ($msg = "$msg\n");
   warn("\n$msg$help\n");

   warn("Scripts in default location ($installDir/R/RExecutor/):\n" 
       .`cd $installDir/R/RExecutor;ls -1` ."\n");

  exit(1);
}


#Write critical command text to $commandLogFile
sub printToCommandLogfile
{  #Function variable declaration and assignment to the local variables.
   my( $commandLogFile, $LogText ) = @_;
   open( COMMANDLogFile, ">>$commandLogFile" ); 
   print( COMMANDLogFile "$LogText \n\n" );
   close( COMMANDLogFile );
}

sub readFile
{  my($rScript) = @_;
   my($txt) = "";
   
   open(IN,$rScript) || die $!;
   while($_=<IN>)
   {  $txt .= $_;
   }
   close(IN);
   return $txt;
}
