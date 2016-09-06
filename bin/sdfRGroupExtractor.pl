#!/usr/bin/env perl
use warnings;

$script=$0;
$script =~ /^\// || ( $script="$ENV{PWD}/$script" );
our($installDir) = $script =~ /(.*)\/[^\/]+$/;

use Getopt::Long;
use File::Copy;
use File::Path;


my $PROGNAME = "sdfRGroupExtractor.pl\n";
my $UPD_DATE = "08/28/2015";
my $HELP_TEXT = <<___help;
$PROGNAME
Fragment the input molecules as specified by the transformations
NOTE: SMIRKS will be applied before SMARTS if -smirksFile and -smartsFile are specified. 
(Last Update: $UPD_DATE)
   -in <arg>...Supported file format: .sdf
   -out <arg>...Supported file format: .sdf
          Fields with RGroup SMILES and ID are added.
   -rgOutFormat <arg>...Supportd RGroup output file format sdf, tab, and csv (comma-delimited)
         If not specified, then it is csv
   -smirksFile <arg>...File with SMIRKS and corresponding name separated by tab.
         Use [U+1], [U+2], etc. to specify the attachment point.
         Names for SMIRKS should not contain spaces
   -smartsFile <arg> ... File with SMARTS and corresponding name separated by tab.
         Use [U+1], [U+2], etc. to specify the attachment point.
         Names for SMARTS should not contain spaces
   -makeH <explicit|implicit>...If not specified, hydrogens in molecules will be made explicit.
   -transform <once|all>...once: used the first applicable transformation.
         once if not specified
   -debug...Print additional info and do not delete temporary files.
   -h...Output the help page
___help


my $USER     = $ENV{'USER'};
my $MYPREFIX = "_RGE__$$";
my $MYID_TAG      = "_internalID_";
my $MYTRANSBY_TAG = "transformedBy";   # Tag as output by sdfTransformer.csh
my $RG_SMILES_TAG = "RG SMILES";
my $RG_TYPE_TAG   = "RG Type";
my $RG_ID_TAG     = "RG ID";

# Set based on user's option specification
my $rgOutFormat = "csv";
my $makeHFlag   = "explicit";
my $TransFlag   = "once";


GetOptions( 'in=s' => \$input, 
            'out=s' => \$output, 
            'smirksFile=s'  => \$smirksFile,
            'smartsFile=s'  => \$smartsFile,
            'rgOutFormat=s' => \$rgOutFormat, 
            'makeH=s' => \$makeHFlag,
            'transform=s' => \$TransFlag,
            'debug' => \$debug,
            'h' => \$help );
if( $help )
{  warn( "$HELP_TEXT");
   exit();
}
warn( "${PROGNAME}(Last Update: $UPD_DATE)\n\n" );


$errMessage = "";
if( !$input )                                                                   
{  $errMessage .= "-in must be specified\n"; 
} elsif( $input !~ /.*\.sdf(\.gz)?/ )
{   $errMessage .= "input type must be .sdf\n";
} else
{  if( $input ne ".sdf" )
   {  -e $input || die "$input does not exists";
   } 
}
if( !$output )                                                                   
{  $errMessage .= "-out must be specified\n";
} elsif( $output !~ /.*\.sdf(\.gz)?/ )
{   $errMessage .= "output type must be .sdf\n";
}
if( $rgOutFormat !~ /sdf|tab|csv/ )
{  $errMessage .= "-rgOutFormat must be sdf, tab, or csv\n";
}
if( !$smirksFile && !$smartsFile )
{  $errMessage .= "smirksFile or smartsFile must be specified\n"; 
}
if( $makeHFlag !~ /explicit|implicit/ )
{  $errMessage .= "-makeH must be explicit or implicit\n";
} elsif( $makeHFlag =~ /implicit/ )
{  $makeHFlag = "-makeHImplicit";
} else
{  $makeHFlag = "-makeHExplicit";
}
if( $TransFlag !~ /once|all/ )
{  $errMessage .= "-transform must be once or all\n";
} elsif( $TransFlag =~ /all/ )
{  $TransFlag = "";
} else
{  $TransFlag = "-transformOnce";
}

$errMessage && die( "ERROR\n----\n$errMessage \n-----\n$HELP_TEXT" );
writeToLogFile( "$rgOutFormat" );


# Add internal record id and apply transformations
# Following files are created
#  ${MYPREFIX}_input.sdf   [Replace existing molfile title with MYID]
#  ${MYPREFIX}_tf.sdf      [Molfile of fragments in place of the original molfile]
#  ${MYPREFIX}_tf.tab      [Mapping of internal id to SMIRKS used to transform]
$readInputCom = "cat $input";
if( $input =~ /.*\.sdf\.gz/ )
{  $readInputCom = "gunzip -c $input";
}elsif( $input =~ /^\.sdf/ )
{  $readInputCom = "cat -";
}
my $transOption = "";
if( $smirksFile )
{  $transOption = "-trans $smirksFile";
}
if( $smartsFile )
{  $transOption = "-scaffold $smartsFile";
}

$com = <<'___COMS';
#readInputCom# \
   | sdfTagTool.csh -in .sdf -out .sdf -addCounter -counterTag #MYID_TAG# \
   | sdfTagTool.csh -in .sdf -out .sdf -title #MYID_TAG# -remove "counter" \
   | tee  #MYPREFIX#_input.sdf \
   | sdfTransformer.csh -in .sdf -out .sdf \
      #makeHFlag# #TransFlag# #transOption# \
   | tee  #MYPREFIX#_tf.sdf \
   | sdf2Tab.csh -in .sdf -tags "#MYID_TAG#|#MYTRANSBY_TAG#" \
   | perl -ne 'chomp; if( /.+\t.+/ ){ print "$_\n"; }' \
   > #MYPREFIX#_tf.tab
___COMS
$com =~ s/#readInputCom#/$readInputCom/g;
$com =~ s/#MYID_TAG#/$MYID_TAG/g;
$com =~ s/#MYPREFIX#/$MYPREFIX/g;
$com =~ s/#makeHFlag#/$makeHFlag/;
$com =~ s/#TransFlag#/$TransFlag/;
$com =~ s/#transOption#/$transOption/;
$com =~ s/#MYTRANSBY_TAG#/$MYTRANSBY_TAG/;
writeToLogFile( "\n\n$com" );
$status = system( $com ) == 0 or die "$_";


#Replace "*+" with "R", generate the RGroup type and id
# Following files are created
#  ${MYPREFIX}_fragAll.tab       [Map RGroup to internal id]
#  ${MYPREFIX}_fragUnique.tab    [Unique list of RGroups, i.e. map RGroup to RGroup id]
#  ${MYPREFIX}_fragTypes.tab     [RGroup types for compiling individual merge files]
#  ${MYPREFIX}_fragMerged.tab    [RGroups id and SMILES mapped to internal id]
$com = <<'___COMS';
perl -e 'print "#RG_SMILES_TAG#\t#MYID_TAG#\t#RG_ID_TAG#\n";' > #MYPREFIX#_fragAll.tab
perl -e 'print "#RG_SMILES_TAG#\t#RG_TYPE_TAG#\t#RG_ID_TAG#\n";' > #MYPREFIX#_fragUnique.tab

cat #MYPREFIX#_tf.sdf \
   | sdfMolSeparator.csh -in .sdf -out .smi \
   | perl -pe 's/\[UH\+(\d*)\]/\[H\]\[U\+$1\]/' \
   | perl -pe 's/(.+) (.+)/$1\t$2/' \
   | perl -pe 's/\[U\+\]/\[U\+1\]/g' \
   | perl -ne '/U\+\d/ && print' \
   | perl -pe 'chomp; @a=/(U\+\d+)/g; @a=sort(@a); $_=$_."\t".join("_",@a)."\n";' \
   | tee -a #MYPREFIX#_fragAll.tab \
   | perl -pe 'chomp; @_=split(/\t/); $_="$_[0]\t$_[2]\n";' \
   | sort | uniq \
   | perl -pe 'chomp; @_=split(/\t/); $c=$_[1]."_".$cnt{$_[1]}++; $c=~s/\+//g; $_=join("\t",@_)."\t".$c."\n";' \
   | tee -a #MYPREFIX#_fragUnique.tab \
   | perl -pe 'chomp; @_=split(/\t/); $_="$_[1]\n";' \
   | sort | uniq \
   > #MYPREFIX#_fragTypes.tab

tabTabMerger.pl -add #MYPREFIX#_fragUnique.tab #MYPREFIX#_fragAll.tab \
   | perl -pe 'chomp; @_=split(/\t/); $_=$_[1]."\t".$_[3]."\t".$_[5]."\t".$_[4]."\n";' \
   > #MYPREFIX#_fragMerged.tab
___COMS
$com =~ s/#RG_SMILES_TAG#/$RG_SMILES_TAG/g;
$com =~ s/#RG_ID_TAG#/$RG_ID_TAG/g;
$com =~ s/#RG_TYPE_TAG#/$RG_TYPE_TAG/g;
$com =~ s/#MYPREFIX#/$MYPREFIX/g;
$com =~ s/#MYID_TAG#/$MYID_TAG/g;
writeToLogFile( "\n\n$com" );
$status = system( $com );
$status and die "$!";


#Compile the files for adding to the input file
# Following files are created
#  ${MYPREFIX}_fragMerged.tab       [RGroups id and SMILES mapped to internal id]
$sepCom = <<'___sepCOMS';
cat #MYPREFIX#_fragMerged.tab \
   | perl -ne 'if( /\t#RG_TYPE#\n|#MYID_TAG#/ ){ print $_;}' \
   | perl -pe 's/RG /#RG_TYPE_TAG# /g' \
   | perl -pe 'chomp; @_=split(/\t/); $_="$_[0]\t$_[1]\t$_[2]\n";' \
   > #MYPREFIX#_#RG_TYPE_TAG#.tab
___sepCOMS
$sepCom =~ s/#MYPREFIX#/$MYPREFIX/g;
$sepCom =~ s/#MYID_TAG#/$MYID_TAG/g;

$addCom = <<___addCOMS;
   | sdfTabMerger.csh -out .sdf -sdf .sdf -tab ${MYPREFIX}_#RG_TYPE_TAG#.tab \\
         -outAll -quiet -addEmptyValues -mergeTag "$MYID_TAG" -mergeCol "$MYID_TAG" \\
___addCOMS

#Compile the files for adding to the input file
$com = "";
$mergeCom = "";
open( FILE, "${MYPREFIX}_fragTypes.tab" ) or die "${MYPREFIX}_fragTypes.tab \n $!";
while( $line = <FILE> )
{
   chomp( $line );
   $rgTypeTag = $line;
   $rgTypeTag =~ s/\+//g; # Names look nicer
   $line =~ s/\+/\\+/g;   # Make RG_TYPE string ready for RegExp match
   $com .= "\n".$sepCom;
   $com =~ s/#RG_TYPE#/$line/g;
   $com =~ s/#RG_TYPE_TAG#/$rgTypeTag/g;
   $mergeCom .= $addCom;
   $mergeCom =~ s/#RG_TYPE#/$line/g;
   $mergeCom =~ s/#RG_TYPE_TAG#/$rgTypeTag/g;
}
close( FILE );

writeToLogFile( "\n\n$com" );
$status = system( $com );
$status and die "$!";

$com = <<___COMS;
cat ${MYPREFIX}_input.sdf \\
   | sdfTabMerger.csh -out .sdf -sdf .sdf -tab ${MYPREFIX}_tf.tab \\
         -outAll -quiet -addEmptyValues -mergeTag "$MYID_TAG" -mergeCol "$MYID_TAG" \\
$mergeCom   | sdfTagTool.csh -in .sdf -out $output -remove "$MYID_TAG"
___COMS
writeToLogFile( "\n\n$com" );
$status = system( $com );
$status and die "$!";


$tabCom = <<'___COMS';
cat #MYPREFIX#_fragUnique.tab \
   | perl -ne 'if( /\t#RG_TYPE#\t|#RG_ID_TAG#/ ){ print $_;}' \
   | perl -pe 'chomp; @_=split(/\t/); if(/RG ID/){$t="TITLE";}else{$t=$_[2];} $_="\"$_[0]\"#SEP#\"$t\"#SEP#\"$_[2]\"#SEP#\"$_[1]\"\n";' \
___COMS
$tabCom =~ s/#MYPREFIX#/$MYPREFIX/g;
$tabCom =~ s/#RG_ID_TAG#/$RG_ID_TAG/g;

$sdfCom = <<'___COMS';
   | perl -pe 'chomp; @_=split(/,/); $_="$_[0],$_[1],$_[2],$_[0],$_[3]\n";' \
   | sdfTagTool.csh -in .csv -out .sdf \
___COMS

my $sep = "\\t";
if( $rgOutFormat eq "csv" )
{  $sep = ",";
} elsif( $rgOutFormat eq "sdf" )
{  $sep = ",";
   $tabCom .= $sdfCom;
}
$tabCom .= "   > LIST_#RG_TYPE_TAG#.#rgOutFormat#";
$tabCom =~ s/#SEP#/$sep/g;
$tabCom =~ s/#rgOutFormat#/$rgOutFormat/g;
writeToLogFile( "\n\n$tabCom" );

open( FILE, "${MYPREFIX}_fragTypes.tab" ) or die "${MYPREFIX}_fragTypes.tab \n $!";
while( $line = <FILE> )
{
   chomp( $line );
   $rgTypeTag = $line;
   $rgTypeTag =~ s/\+//g; # Names look nicer
   $line =~ s/\+/\\+/g;   # Make RG_TYPE string ready for RegExp match
   $com = $tabCom;
   $com =~ s/#RG_TYPE#/$line/g;
   $com =~ s/#RG_TYPE_TAG#/$rgTypeTag/g;
   writeToLogFile( "\n\n$com" );
   $status = system( $com );
   $status and die "$!";
}
close( FILE );


#Remove tmp files
if( !$debug )
{  $com = "\\rm ${MYPREFIX}_*";
   writeToLogFile( "\n\n$com" );
   system( $com ) == 0 or die "$!";
}




#--------------------------------- Utility Methods ---------------------------------#
sub writeToLogFile
{  my( $text ) = @_;
   open( LOGFILE, ">>sdfRGExtractor_Log.txt" );
   print( LOGFILE "${text}\n\n" );
   close( LOGFILE );
}

