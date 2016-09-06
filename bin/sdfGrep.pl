#!/usr/bin/env perl
use warnings;
$script=$0;
$script =~ /^\// || ( $script="$ENV{PWD}/$script" );
our($installDir) = $script =~ /(.*)\/[^\/]+$/;

use File::Copy;
use Getopt::Long;

my($omegaOpt ) = " -ewindow 20 -maxconfs 500 ";
my($szybkiOpt) = " -grad_conv .001 -max_iter 2000 ";

my( @constrVals ) = ("0.2", "0.6", "1.0", "1.4", "1.8");

$use = <<USE;

sdfGrep.pl [-v][-i][-s][-m] regEx sdf....
  -v .............. return non matching records
  -i .............. ignore case
  -m .............. ^ and \$ match line start/end
  -s .............. . matches newlines
The regEx will be applied to the each full sdf record incuding
molblock, fieldnames and field values.
For deatials on regex type "man perlre".

USE

my( $invert, $ignoreC, $singleLine, $multiLine) = ("", "", "", "");
GetOptions("v" =>\$invert,
           "i" => \$ignoreC,
           "m" => \$singleLine,
           "s" => \$multiLine,
          ) || die $use;
$#ARGV == -1 && die "Missing pattern\n$use";

my($re) = shift;
my($modifier) = "";
$ignoreC &&    ($modifier = "i");
$singleLine && ($modifier .= $singleLine );
$multiLine &&  ($modifier .= $multiLine );

$modifier && ($re = "(?$modifier)$re");

while($_ = <>)
{  $rec .= $_;
   if( /\$\$\$\$/ )
   {  if( (!$invert && $rec =~ /$re/) || ($invert && $rec !~ /$re/) )
      {  print $rec;
      }

      $rec = "";
   }
}
