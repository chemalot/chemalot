#!/bin/csh -f
#

set main=autocorrelator.apps.SdfSdfMerger

set script=$0
if( "$script" !~ "/*" ) set script=$PWD/$script
set installDir=$script:h

source $installDir/starter_csh

