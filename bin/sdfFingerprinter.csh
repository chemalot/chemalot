#!/bin/csh -f
#

set main=com.aestel.chemistry.openEye.fp.Fingerprinter
if($?JAVAXMX) then
  set XMX=$JAVAXMX
else
   set XMX=30G
endif


set script=$0
if( "$script" !~ "/*" ) set script=$PWD/$script
set installDir=$script:h

source $installDir/starter_csh

