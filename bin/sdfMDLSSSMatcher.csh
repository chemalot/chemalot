#!/bin/csh -f
#

set main="com.genentech.chemistry.openEye.apps.SDFMDLSSSMatcher"
if($?JAVAXMX) then
  set XMX=$JAVAXMX
else
   set XMX=2G
endif


set script=$0
if( "$script" !~ "/*" ) set script=$PWD/$script
set installDir=$script:h

source $installDir/starter_csh

