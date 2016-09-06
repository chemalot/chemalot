#!/usr/bin/env perl
use warnings;
$script=$0;
$script =~ /^\// || ( $script="$ENV{PWD}/$script" );
our($installDir) = $script =~ /(.*)\/[^\/]+$/;

use Getopt::Long;
use Cwd;

$use = "\n"
  ."tabTabMerger.pl -add addFile.tab [inFile.tab] [>outFile.tab]\n"
  ."Will read read addFile into a hashmap and join it to inFile\n"
  ."The join will be based on the first column in both files\n"
  ."The first row in both files is assumed to contain headers\n"
  ."Records in the inFile which have no record in the add file will\n"
  ."  be returned unchanged. Records in the addFile which have no\n"
  ."  record in the inFile will be silently ignored.\n"
  ."Files MUST have unix linefeeds!\n"
  ."\n";

my( $addFile ) = ('');
GetOptions("-addFile=s" => \$addFile );

$addFile || die "$use\n";
$#ARGV > -1 && !-r $ARGV[0] && die "$use\n";

open( ADD, "$addFile" ) || die "Could not open $addFile $!\n";
$header = <ADD>;
chomp($header);
@header = split( /\t/, $header );
while($_=<ADD>)
{  chomp;
   @_ = split(/\t/);
   $key = $_[0];
   if( $#_ > $#header )
   {   die "line contains more values than header:\n$_\n";
   }

   if( ! $addMap{$key} )
   {  $addMap{$key} = $_;
   } else
   {  warn "repeat found for key $key, using first.\n";
   }
}
close(ADD);


$_ = <>;
chomp;
@inHeader = split(/\t/);
print( "$_\t$header\n" );
while(<>)
{  chomp;
   if( ! $_ )
   {  print "\n";   # empty line
      next;
   }
   @_ = split(/\t/,$_,-1);

   if( $#_ > $#inHeader )
   {  die "Row has more values than header:\n$_\n";
   }else
   {  $_ = $_ . ("\t" x ($#inHeader - $#_));
   }

   $key = $_[0];

   if( defined( $addMap{$key} ) )
   {  $add = $addMap{$key};
   }else
   {  $add = "";
   }

   print "$_\t$add\n";
}
   
