#!/bin/tcsh -f

if("$argv" == '-h') then
   echo 'call using similar to one of: '
   echo '   set n=XXX.g; set nCPU=2; bsub -q medium -J $n:r -oo $n:r%J.qlog -n $nCPU -R "rusage[mem=5] span[hosts=1]" gaussian.csh $n'
   echo '   set n=XXX.g; set nCPU=4; bsub -q medium -J $n:r -oo $n:r%J.qlog -n $nCPU -R "rusage[mem=3] span[hosts=1]" gaussian.csh $n'
endif

if( -e /gne/research/scratch/$user ) then
   setenv GAUSS_SCRDIR /gne/research/scratch/$user
else
   setenv GAUSS_SCRDIR /gne/research/scratch/smdi
endif
if( ! -e $GAUSS_SCRDIR/gaussian ) mkdir $GAUSS_SCRDIR/gaussian
setenv GAUSS_SCRDIR $GAUSS_SCRDIR/gaussian/$HOSTNAME
if( ! -e $GAUSS_SCRDIR ) mkdir $GAUSS_SCRDIR

(find $GAUSS_SCRDIR -mtime +30 -type f -not -iname =readme\* -delete; \
 find $GAUSS_SCRDIR -mindepth 1 -mtime +2  -type d -empty -delete) \
|& perl -ne 'warn $_'

limit coredumpsize 0

if( $#argv == 1 ) then
   set inFile="$1"
   g09 <"$inFile" >> "$inFile:r.out"
else
   g09
endif
