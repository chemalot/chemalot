#!/usr/bin/env perl
use warnings;
use Getopt::Long;
$use = "g2sdf.pl [-mem nnGB] [-FNTag tagName|-FNPrefix prefix] [-template gausianTemplate.g]\n"
      ."         [-fixTorAtomTag tag] sdfFiles\n"
      ."Writes one file per sdf record named by the value of tagname.\n"
      ."If -FNPrefix is given then the output .g files are numbered.\n"
      ."fixTorAtomTag . Specify name for sdf field containing , separated list of fixed torsion atoms (0 based)\n"
      ."If neither FNTag nor FNPrefix is given then all out is writen to inName.g\n"
      ."   If input is from stdein (no inName) output will be to stdout\n"
      ."If tagname is TITLE the TILE is used.\n"
      ."Default template is AM1 optimization.\n"
      ."Tamplate must contain #XYZ# string\n"
      ."The default memis 10GB\n";
my($script)=$0;
$script =~ /^\// || ( $script="$ENV{PWD}/$script" );
my($installDir) = $script =~ /(.*)\/[^\/]+$/;
my($guasTemplateDir) = "$installDir/data/gaussian";

my(  $tFile, $ofnTag,$fnPrefix, $mem, $fixTorAtomTag ) = ('', '', '', '10GB', '');
GetOptions("FNTag=s" =>\$ofnTag,"FNPrefix=s"=>\$fnPrefix, "template=s" =>\$tFile,
           "mem=s" =>\$mem, "fixTorAtomTag=s" => \$fixTorAtomTag)
  || &exitHelp();

#read template
my($template) = '';
$tFile || ($tFile = "$guasTemplateDir/default.g");
if( -e $tFile )
{  $template = `cat $tFile`;
}else
{  if( -e "$guasTemplateDir/$tFile" )
   {  $template = `cat "$guasTemplateDir/$tFile"`;
   }else
   {  die "tempalte $template not found!\n";
   }
}

if( $fixTorAtomTag )
{  if( $template !~ /#FIX#/ )
   {  die "Template needs to conatin #FIX# placeholder for fix atom specification\n";
   }
   if( $template !~ /ModRedundant/i )
   {  die "Template must use ModRedundant option for gaussian to use fixed atoms\n";
   }
}

# loop over input files
my($record,$fixTorAtoms) = ('','');
my( $num, $lastFName ) = ( 1, '');
while(<>)
{  if( $lastFName ne $ARGV )
   {  if( $record && $lastFName )
      {  # last file did not end in $$$$
         $fName = &cleanFName($lastFName);
         &outSDF( $fName, $num++, $record );
         $record = '';
         $fName = '';
         $fixTorAtoms = '';
      }
       $lastFName =  $ARGV;
   }

   $record .= $_;
   if( ! $fName && $ofnTag eq 'TITLE' )
   {  $fName = $_;
   }

   if( m/<$ofnTag>/ )
   {  $fName = <>;

   } elsif( m/<$fixTorAtomTag>/ )
   {  $fixTorAtoms = <>;
      $fixTorAtoms = &convertToAtList($fixTorAtoms);
      $fixTorAtoms = "D $fixTorAtoms F";

   } elsif( $_ =~ m/^\$\$\$\$/ )
   {  if( !$ofnTag && ! $fnPrefix ) { $fName = &cleanFName($ARGV); };
      &outSDF(  $fName, $num++, $record );
      $record = '';
      $fName = '';
   }
}
if( $record )
{  if( !$ofnTag && ! $fnPrefix ) { $fName = &cleanFName($ARGV); };
   &outSDF(  $fName, $num++, $record ); 
}

# remove trailing .sdf
sub cleanFName
{  my( $fName ) = @_;

   $fName =~ s/\.[^.]+$//;
   return $fName;
}


sub outSDF
{  my( $fName, $num, $record ) = @_;
   my( $gCom, $appendMode, $xyz, @xyz ) = ( '', '>', '' );

   if( $fName ) 
   {  $fName =~ s/\s$//g; 
      $fName =~ s/\s/_/g;
      $fName =~ s/\.sdf$//ig;
   }

   if( $fnPrefix )
   {  $fName = "$fnPrefix$num";
      $appendMode = '>';
   }


   $xyz = `$ENV{OBABEL} -isdf -ogau <<sdf2g\n$record\nsdf2g`;
   @xyz = split(/\n/, $xyz);
   $xyz = join("\n", @xyz[4 .. $#xyz]);
   $gCom = $template;
   $gCom =~ s/#XYZ#/$xyz/g;
   $gCom =~ s/#FName#/$fName/g;
   $gCom =~ s/#mem#/$mem/g;
   $gCom =~ s/#FIX#/$fixTorAtoms/g;

   if( $fName eq "-" && !$ofnTag && ! $fnPrefix )
   {  # input was from stdin ("-") and no ofnTag or FNPrefix was given
      print $gCom;
   } else
   {  if( $fName )
      {  $fName =~ s/\s$//g;
         $fName =~ s/\s/_/g;
         $fName =~ s/\.sdf$//ig;
      }
   
      if( $fnPrefix )
      {  $fName = "$fnPrefix$num";
         $appendMode = '>';
      }
      open(OFILE, "$appendMode$fName.g") || die "$appendMode$fName.g: $!";
      print(OFILE $gCom);
      close(OFILE);
   }
}

sub exitHelp
{  warn("$use\nAvailable default templates (check: $guasTemplateDir):\n");
   while( <$guasTemplateDir/*.g>)
   {  s/.*\///;
      warn("$_\n");
   }
   die "\n";
}

sub convertToAtList
{  my(@atList) = split(/ /,$_[0]);
   my($at);
   for $at (@atList) { $at++; }
   return join(" ", @atList);
}

