#!/bin/csh -f

set help=0
if($#argv == 0) then
   set help=1
else 
   if( $argv[1] == '-h') set help=1
endif
if( $help == 1) then
  cat <<COMS|perl -ne 'warn($_)'
  sdfStructureTagger.csh -in .sdf -out .sdf
     Will add counts of the 166 Bit MACCS keys to the input file.
COMS
  exit 1
endif

sdfStructureTagger.csh -smarts $AESTEL_DIR/config/fp/MACCSKeys.txt -sets MACCS \
                       -output_exists $*:q
