#!/bin/csh -f
#

set main="com/genentech/chemistry/openEye/apps/SDFTopologicalIndexer"

set script=$0
if( "$script" !~ "/*" ) set script=$PWD/$script
set installDir=$script:h
if( ! $?javaOpts ) then
   #optimized for large object creation
   set javaOpts="-XX:NewRatio=4 -Xmx2G"
endif

source $installDir/starter_csh

