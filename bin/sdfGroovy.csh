#!/bin/csh -f
#

set main=autocorrelator.apps.SDFGroovy
if($?JAVAXMX) then
  set XMX=$JAVAXMX
else
   set XMX=30G
endif


set script=$0
if( "$script" !~ "/*" ) set script=$PWD/$script
set installDir=$script:h
if( $#argv > 0 && "$1" !~ "-*" ) then
   set argv=(-f $argv:q)
endif

source $installDir/starter_csh

