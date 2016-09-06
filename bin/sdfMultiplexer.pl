#!/usr/bin/env perl
use warnings;


$class="com.genentech.chemistry.tool.sdfMultiplexer.SDFMultiplexer";

$script=$0;
$script =~ /^\// || ( $script="$ENV{PWD}/$script" );
our($installDir) = $script =~ /(.*)\/[^\/]+$/;
our($cp) = "$installDir/../lib/*";

@javaOpts=("-Xmx5G", "-Xms256m");
if( defined $ENV{"OE_LIB_DIR"} ) 
{  push(@javaOpts, "-Doejava.libs.path=$ENV{OE_LIB_DIR}"); 
}

warn "$ENV{HOST}\n";
@cmdList = ("java",
  @javaOpts, "-cp", $cp, $class, @ARGV );
#warn join(" ", @cmdList) . "\n";
exec( @cmdList );
