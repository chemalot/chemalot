#!/bin/csh -f

#author JW Feng
set main="com/genentech/application/calcProps/SDFCalcProps"

set script=$0
if( "$script" !~ "/*" ) set script=$PWD/$script
set installDir=$script:h
set basename = $script:t:r
set libDir=$installDir/../lib

if( ! $?XMX )      set XMX=1G
if( ! $?javaOpts ) set javaOpts=-Xmx$XMX
set javaOpts="$javaOpts -XX:+UseMembar"

set jCom=(java $javaOpts -cp "'$libDir/*'" $main $*:q)
set command="`$jCom:q`"
if( $status != 0 ) exit 1

#execute the command
#perl -e "warn '$command:q'"
eval $command:q
