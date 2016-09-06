#!/usr/bin/env perl
use warnings;
$script=$0;
$script =~ /^\// || ( $script="$ENV{PWD}/$script" );
our($installDir) = $script =~ /(.*)\/[^\/]+$/;


use Getopt::Long;
$use = "sdfSliceByRe.pl (-after re |-v -filter re) <inFile\n"
  ."   -after re .. print all records following the first record matching re\n"
  ."   -filter re . print any record matching re\n"
  ."   -v ......... print any records not matching '-filter re'\n"
  ."   Piping seems not to work on Windows\n";


my( $afterRE ) = '';
my( $filterRE ) = '';
my( $invert ) = 0;
GetOptions("after=s" => \$afterRE, "filter=s" => \$filterRE, "v" => \$invert );
if( $afterRE )
{  &afterRE( $afterRE );
   exit(0);
}

if( $filterRE )
{  &filterRE( $filterRE, $invert );
   exit(0);
}

die( $use );

sub filterRE
{  my( $filterRE, $invert ) = @_;
   my( $mol, $matched );

   while($_ = <STDIN>)
   {  $mol .= "$_";
      if( /^\$\$\$\$/ )
      {  if( ($matched && ! $invert) || (!$matched && $invert) )
         {  print $mol;
         }
         $matched = 0;
         $mol = '';
      }else
      {  if( /$filterRE/ )
         {  $matched = 1;
         }
      }
   }

   exit(0);
}
   

sub afterRE
{  my($afterRE) = @_;
   $afterRE || die $use;
   $#ARGV > 0 && die $use;

   $file = shift(@ARGV);
   $file && (open(STDIN, $file) || die "$!");

   while($_ = <STDIN>)
   {  /$afterRE/ && last;
   }

   while($_ = <STDIN>)
   {  /^\$\$\$\$/ && last;
   }

   while($_ = <STDIN>)
   {  print;
   }
   exit(0);
}

