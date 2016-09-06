#!/usr/bin/env perl
use warnings;

use Getopt::Long;

# help not taken from java code because perl does handle directory lookup better
$use = "tabExport.pl sqlFileName sqlTaskName [-noHeader] [-o outName] [-newLineReplacement] params\n"
      ."  sqlFileName may be relative to \$AESTEL_DIR/config/exporter\n"
      ."  -newLineReplacement <arg>   If given newlines in fields will be replaced\n";
$class="com.genentech.retrival.tabExport.TABExporter";
$script=$0;
$script =~ /^\// || ( $script="$ENV{PWD}/$script" );
our($installDir) = $script =~ /(.*)\/[^\/]+$/;
our($cp) = "$installDir/../lib/*";

$#ARGV >= 1 || die "$use\n";

$sqlFile = shift;
$sqlName  = shift;

-e $sqlFile || ( $sqlFile = "$ENV{AESTEL_DIR}/config/exporter/$sqlFile" );
(-e $sqlFile && ! -d $sqlFile) || ( $sqlFile = "$sqlFile/sql.xml" );

$script=$0;
$script =~ /^\// || ( $script="$ENV{PWD}/$script" );

@javaOpts=();
if( defined $ENV{"OE_LIB_DIR"} ) { push(@javaOpts, "-Doejava.libs.path=$ENV{OE_LIB_DIR}"); }

warn "$ENV{HOST}\n";
@cmdList = ("java",
  "-Djava.security.egd=file:///dev/urandom",
  @javaOpts, "-cp", $cp, $class, 
  "-sqlFile", $sqlFile, "-sqlName", $sqlName, @ARGV );
warn join(" ", @cmdList) . "\n";
exec( @cmdList );
