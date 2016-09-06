#!/bin/csh -f
#

set main="com.genentech.chemistry.openEye.apps.SDFMolSeparator"
if($?JAVAXMX) then
  set XMX=$JAVAXMX
else
   set XMX=1G
endif


set script=$0
if( "$script" !~ "/*" ) set script=$PWD/$script
set installDir=$script:h

source $installDir/starter_csh

