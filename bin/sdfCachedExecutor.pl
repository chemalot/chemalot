#!/usr/bin/env perl
use warnings;
$script=$0;
$script =~ /^\// || ( $script="$ENV{PWD}/$script" );
our($installDir) = $script =~ /(.*)\/[^\/]+$/;

use File::Copy;
use Getopt::Long;
use File::Temp qw/ tempfile /;

$use = <<USE;
sdfCachedExecutor.pl -in .sdf -out .sdf [-inId id] [-cacheDir dir] -cmd cmd
         [-ignorePat regex] [-shell sh]
  -inId ....... external identifier for input file used on top of checksum
                if checksum, id and cmd are identical the chaced resutls are used
  -cacheDir ... directory to store cache files def: /tmp/sdfCache/
  -ignorePat .. ignore lines matching pattern for checksum def: '^  -OEChem-\\d+'
  -shell ...... Overwrite /bin/sh as the executing shell
  -cmd ........ command to execute, must work as pipe.
USE


our($tmpBase) = "/tmp/sdfCache.$$";
our($debug) = "";
mkdir($tmpBase);
my( $com, $checkSum, $cachedFile, $fType) = ("", "");
my( $in, $out, $cacheDir, $ignorePat, $inId, $cmd, $shell) 
  = ("", "", "", "", "", "");
GetOptions("in=s" =>\$in, 
           "out=s" =>\$out,
           "inId=s" =>\$inId,
           "cmd=s" =>\$cmd,
           "cacheDir=s" =>\$cacheDir,
           "ignorePat=s" => \$ignorePat,
           "shell=s" => \$shell,
           "debug" => \$debug
           ) || die $use;
$#ARGV == -1 || die $use;

$in && $out || die "in and out needed\n\n$use\n";
$cmd || die "cmd missing\n\n$use\n";

($fType) = $in =~ /(\.\w+)$/;
$fType || die "input file has no extention\n\n$use\n";
$fType ne "CSum" || die "input file type may not be CSum\n\n$use\n";

$ignorePat || ($ignorePat = "^  -OEChem-\\d+");
$cacheDir || ($cacheDir = "/tmp/sdfCache");
mkdir($cacheDir);
$inId     || ($inId = "noInId");

($in, $checkSum) = &getCheckSum($in, $cmd, $inId, $ignorePat);
$cachedFile = &getCachedFile($cacheDir, $checkSum, $cmd, $inId, $fType);

if( $cachedFile )
{  warn( "Using cached data: $cachedFile.gz\n" );
   if( $out =~ /^\.\w+$/ )
   {  system("gunzip -c $cachedFile.gz");
   } else
   {  system("gunzip -c $cachedFile.gz >$out");
   }

   if( !$debug )
   {  unlink(<$tmpBase/*>);
      rmdir($tmpBase);
   }
   exit(0);
}

my($cacheFileBase) = &createCacheEntry($cacheDir, $checkSum, $cmd, $inId, $fType);

$com = "cat '$in'|$cmd|tee $cacheFileBase$fType";
if( $out !~ /^\.\w+$/ )
{  $com .= ">$out";
}

# redirect stderr to log file before execution
open(orgSTDERR,">&STDERR") ||die "Could not redirect sdterr: $!";
open(STDERR, ">$cacheFileBase.log") || die "Could not open error file: $cacheFileBase.log: $!";

if( $debug ) { print orgSTDERR "$com\n"; }
# execute 
my(@com) = ($com);
$shell && (@com = ($shell,"-c", $com));
my($status) = system(@com);
$status += system("gzip $cacheFileBase$fType");

open(STDERR, ">/dev/null"); # close(STDERR) withour warning
*STDERR=orgSTDERR;
*orgSTDERR=STDERR;
open(ERROUT, "<$cacheFileBase.log") || die "Could not open error file for reading: $cacheFileBase.log: $!";
while($_=<ERROUT>) { warn $_; }
close(ERROUT);

if( !$debug )
{  unlink(<$tmpBase/*>);
   rmdir($tmpBase);
}
if($status != 0 )
{  unlink("$cacheFileBase.CSum");
   unlink("$cacheFileBase$fType");
   unlink("$cacheFileBase$fType.gz");
   exit($status);
}

# no make cache valid by adding checksum
open(CACH, ">>$cacheFileBase.CSum");
print CACH "$checkSum\n";
close(CACH);




# create *.CSum file with checkSum, return filename for cached data
sub createCacheEntry
{  my($cacheDir, $checkSum, $cmd, $inId) = @_;
   my( $shortSum, $cacheFile);

   $cmd =~ s/[\n\r\t\a]/ /g;
   $shortSum = substr($checkSum,0,35);
   ($fh, $cacheFile) = tempfile( "$shortSum.XXXXX", SUFFIX => ".CSum",
                        DIR => $cacheDir, unlink_on_destroy => 0, UNLINK => 0);
   print $fh "$cmd\n";
   print $fh "$inId\n";
   close($fh);

   # tempfile will make it not group accessible
   my( $accessMod ) = ((~umask) & 0666);
   chmod( $accessMod, $cacheFile );

   $cacheFile =~ s/\.CSum$//;
   return $cacheFile;
}

# check if cached file is available
sub getCachedFile 
{  my( $cacheDir, $checkSum, $cmd, $inId, $fType) = @_;
   my( $shortSum, $baseName, $cacheFile);

   $cmd =~ s/[\n\r\t\a]/ /g;
   $shortSum = substr($checkSum,0,35);
   while($cacheFile = <$cacheDir/$shortSum.*.CSum>)
   {  ($baseName) = $cacheFile =~ /(.+)\.CSum$/;
      $cacheFile = "$baseName$fType";
      if( ! -e "$cacheFile.gz" )
      {  warn "Datafile does not exist: $cacheFile.gz. deleting...\n";
         unlink( "$baseName.CSum" );
         next;
      }

      open(INFFILE, "$baseName.CSum") || die "Could not read $_:$!";
      $_=<INFFILE>; chomp;
      if( ! $cmd eq $_ )      { close(INFFILE); next; }

      $_=<INFFILE>; chomp;
      if( ! $inId eq $_ )     { close(INFFILE); next; }

      $_=<INFFILE>; chomp;
      if( ! $_ )
      {  warn "No Checksum in $cacheFile! deleteing...\n";
         unlink( "$cacheFile.gz" );
         unlink( "$baseName.CSum" );
         next;
      }elsif( ! $checkSum eq $_ ) 
      {  close(INFFILE);
         next;
      }
      
      close(INFFILE);
      return($cacheFile);
   }
   return "";
}


# stream input and compute checksum
# return filename with original input data and chechSum
#        original input data might be in tmp dir if this was streaming data
sub getCheckSum
{  my($in, $cmd, $inId, $ignorePat) = @_;
   my($com, $outInFile, $ignoreCOM, $checkSum, $saveCOM);

   $outInFile = $in;
   $saveCOM = "";
   if( $in =~ /^\.\w+$/ )
   {  $outInFile = "$tmpBase/infile$in";
      $in = "";
      $saveCOM = "|tee $outInFile";
   }

   $ignoreCOM = "";
   $ignorePat && ( $ignoreCOM = "|perl -ne '/$ignorePat/||print' $in" );
   $cmd  =~ s/'/'\\\''/g;   # quoteexisting single qoutes
   $inId =~ s/'/'\\\''/g;   # quoteexisting single qoutes
   $com = <<COMS;
   ( cat $in $saveCOM $ignoreCOM; echo '$inId'; echo '$cmd') |sha512sum
COMS
   if( $debug ) { warn( "$com\n" ); }

  $checkSum = `$com`;
  $? == 0 || die "$!\n";
  $checkSum =~ s/\s.*//;  #remove filename "-" for stdin
  return ($outInFile, $checkSum);
}
