#!/usr/bin/env perl
use warnings;

use Getopt::Long;

$class="com.genentech.retrival.SDFExport.SDFExporter";
$script=$0;
$script =~ /^\// || ( $script="$ENV{PWD}/$script" );
our($installDir) = $script =~ /(.*)\/[^\/]+$/;
our($cp) = "$installDir/../lib/*";


if( $#ARGV >= 1)
{  $sqlFile = shift;
   $sqlName  = shift;

   -e $sqlFile || ( $sqlFile = "$ENV{AESTEL_DIR}/config/exporter/$sqlFile" );
   (-e $sqlFile && ! -d $sqlFile) || ( $sqlFile = "$sqlFile/sql.xml" );
}

@javaOpts=();
if( defined $ENV{"OE_LIB_DIR"} ) { push(@javaOpts, "-Doejava.libs.path=$ENV{OE_LIB_DIR}"); }

warn "$ENV{HOST}\n";
@cmdList = ("java",
  "-Djava.security.egd=file:///dev/urandom",
  @javaOpts, "-cp", $cp, $class, 
  "-sqlFile", $sqlFile, "-sqlName", $sqlName, @ARGV );
#warn join(" ", @cmdList) . "\n";
exec( @cmdList );
