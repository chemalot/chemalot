#!/bin/csh -f

set main="com/genentech/chemistry/openEye/apps/SDFSubRMSD"

set script=$0
if( $script == $script:h ) set script=$PWD/$script
set installDir=$script:h

source $installDir/starter_csh

