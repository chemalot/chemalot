#!/usr/bin/env perl
use warnings;


# Author: Man-Ling Lee  (started on 04/27/2015)
# Copyright 2015 by Genentech

$script=$0;
$script =~ /^\// || ( $script="$ENV{PWD}/$script" );
our($installDir) = $script =~ /(.*)\/[^\/]+$/;

use Getopt::Long;
use File::Copy;
use File::Path;


my $PROGNAME = "sdfRGroupCalcProps.pl";
my $UPD_DATE = "07/18/2015";
my $HELP_TEXT = <<___help;
$PROGNAME
Replace the attachment points with CH3 before properties calculation.
(Last Update: $UPD_DATE)

usage: $PROGNAME [options] <list of space separated properties>
options:
   -in <arg>...input file (any OE filetype),  for sdtin use .type
   -out <arg>...output file (any OE filetype), for stdout use .type
          Fields with RGroup SMILES and ID are added.
   -propColPrefix <arg>...If not the column names will be as output by sdfCalcProps
   -smiCol <arg>...Name of field containing the R-Group SMILES. This SMILES will be
          used for properties calculation.
   -showProps...Show the list of available properties
   -debug...Print additional info and do not delete temporary files.
   -h...Output help page
___help
writeToLogFile( "$PROGNAME ".join( " ", @ARGV) );


# Set based on user's option specification
my $propColPrefix = "";
my $smiCol = "";                               
my $showProps = 0;

GetOptions( 'in=s' => \$input, 
            'out=s' => \$output, 
            'propColPrefix=s' => \$propColPrefix, 
            'smiCol=s' => \$smiCol, 
            'showProps' => \$showProps, 
            'debug' => \$debug,
            'h' => \$help );
my $propCmd = "";
my $rGroupPropsCmd = "";
for $prop (@ARGV)
{  if( $prop =~ /^rg_/)
   {  $rGroupPropsCmd .= " $prop";
   } else
   {  $propCmd .= " $prop";
   }
}
if( $propCmd )
{  $propCmd = "| sdfCalcProps.csh -in .csv -out .csv $propCmd "
}
if( $rGroupPropsCmd )
{  $rGroupPropsCmd = "| sdfCalcProps.csh -dontFilter -in .csv -out .csv $rGroupPropsCmd "
}


if( $help )
{  warn( "$HELP_TEXT");
   exit();
}
if( $showProps )
{  my $list = `sdfCalcProps.csh |& grep ":" | egrep -v "OE formats|usage:|smdi|resapps|Implemented|aestel"`;
   warn( "$list \n" );
   exit();
}


$errMessage = "";
if( !$input )
{  $errMessage .= "-in must be specified\n";
}
if( !$output )
{  $errMessage .= "-out must be specified\n";
}
if( !$smiCol )
{  $errMessage .= "-smiCol must be specified\n";
}
$errMessage && die( "$errMessage \n$HELP_TEXT" ); 


my $MY_PREFIX = "_RGCP_$$";   #Prefix for the temporary files

#Transformation related variables
my $TRX_METHYL = "[U+][*:1]>>C[*:1]";
my $TRX_H_METHYL = "[UH+]>>CC";        #OpenEyes crashed with -MakeHExplicit on [UH+]


my $inFileType = $input;                                 
$inFileType =~ s/^.*\.([a-z]+)$/$1/;
writeToLogFile( "inFileType ___${inFileType}" );

my $outFileType = $output;
$outFileType =~ s/^.*\.([a-z]+)$/$1/;
writeToLogFile( "outFileType ___${outFileType}" );


if( $propColPrefix )
{  $propColPrefix .= " ";
}
writeToLogFile( "propColPrefix ___${propColPrefix}\n" );


# Setting for case $smiCol is not specified
my $inRGFile = $input;
my $inRGType = $inFileType;
my $tabFile  = $MY_PREFIX."_output.tab";


my $prepInCom = "";
#TODO: Merge this section with the next
#      because now that $smiCol is required both part are always needed
#      However, leave it for time if we want to support when $smiCol is option
#      because the RGroup structure represented by the molfile
if( $smiCol )   # not a list of a unique R-group
{  
   $inRGFile = ".csv";
   $inRGType = "csv";
   
   $prepInCom = <<'___COMS';
sdfTagTool.csh -in #IN_FILE# -out .#IN_TYPE# \
   | tee #MY_PREFIX#_inData.#IN_TYPE# \
   | sdfTagTool.csh -in .#IN_TYPE# -out .csv -rmRepeatTag "#SMI_COL#=2" -reorder "#SMI_COL#" \
   | sdfTagTool.csh -in .csv -out .csv -keep "#SMI_COL#" \
   | perl -pe 'chomp; ($m,$t,$rgs)=split(/,(?=(?:[^"]*"[^"]*")*[^"]*$)/,$_,3); $_="$rgs,$t,$rgs\n";' \
   |
___COMS
   $prepInCom =~ s/#MY_PREFIX#/$MY_PREFIX/g;
   $prepInCom =~ s/#IN_FILE#/$input/g;
   $prepInCom =~ s/#IN_TYPE#/$inFileType/g;
   $prepInCom =~ s/#SMI_COL#/$smiCol/g;
   chomp( $prepInCom );
}


#Calculate properties
# #MY_PREFIX#_input.#IN_TYPE# is generated for debugging purpose only
$com = <<'___COMS';
sdfTagTool.csh -in #IN_FILE# -out .#IN_TYPE# -rename "#SMI_COL#=#MY_PREFIX#_SMILES" \
   | tee #MY_PREFIX#_input.#IN_TYPE# \
   | sdfTagTool.csh -in .#IN_TYPE# -out .csv -keep "#MY_PREFIX#_SMILES" \
   | perl -pe '@_=split(/,(?=(?:[^"]*"[^"]*")*[^"]*$)/,$_,2); $_[0] =~ s/\[(U(H\d?)*)\+\d+\]/\[$1\+\]/g; $_=join(",",@_)' \
   #RGROUP_COMMAND# \
   | sdfTransformer.csh -in .csv -out .csv -makeHImplicit -trans '#TRX_METHYL#' \
   | sdfTransformer.csh -in .csv -out .csv -makeHImplicit -trans '#TRX_H_METHYL#' \
   #PROPS_COMMAND# \
   | sdf2Tab.csh -in .csv \
   | perl -ne '$o=$_; if(/#MY_PREFIX#_SMILES/){ @fs=split(/\t/); $o=join("\t#COL_PREFIX#",@fs);} print $o;' \
   > #TAB_FILE#
___COMS
$com =~ s/#MY_PREFIX#/$MY_PREFIX/g;
$com =~ s/#IN_FILE#/$inRGFile/g;
$com =~ s/#IN_TYPE#/$inRGType/g;
$com =~ s/#RGROUP_COMMAND#/$rGroupPropsCmd/g;
$com =~ s/#TRX_METHYL#/$TRX_METHYL/g;
$com =~ s/#TRX_H_METHYL#/$TRX_H_METHYL/g;
$com =~ s/#PROPS_COMMAND#/$propCmd/g;
$com =~ s/#TAB_FILE#/$tabFile/g;
$com =~ s/#COL_PREFIX#/$propColPrefix/;
$com =~ s/#SMI_COL#/$smiCol/g;

$com = "$prepInCom $com";
writeToLogFile( "$com" );
$status = system( $com );
$status and die "$!";


#Add properties to the input file
$com = <<'___COMS';
sdfTabMerger.csh -sdf  #MY_PREFIX#_inData.#IN_TYPE# -tab #TAB_FILE# \
      -mergeTag "#SMI_COL#" -mergeCol "#MY_PREFIX#_SMILES" \
      -out #OUT_FILE# -outAll -quiet -addEmptyValues
___COMS
$com =~ s/#TAB_FILE#/$tabFile/g;
$com =~ s/#MY_PREFIX#/$MY_PREFIX/g;
$com =~ s/#IN_TYPE#/$inFileType/g;
$com =~ s/#OUT_FILE#/$output/g;
$com =~ s/#SMI_COL#/$smiCol/g;

writeToLogFile( "$com" );
$status = system( $com );
$status and die "$!";

#Remove tmp files
if( !$debug )
{  $com = "rm ${MY_PREFIX}_*";
   writeToLogFile( "\n\n$com" ); 
   system( $com ) == 0 or die "$!";
}



#--------------------------------- Utility Methods ---------------------------------#
sub writeToLogFile
{  my( $text ) = @_;
   open( LOGFILE, ">>sdfRGCalcProps_Log.txt" );
   print( LOGFILE "${text}\n\n" );
   close( LOGFILE );
}





