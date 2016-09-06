#!/bin/csh -f

set pid=$$

echo this tests if the sdfMultiplexer can achieve a speedup
echo by distributing over 10 CPU comapred to runnig on a single CPU
echo in extreem situations this might not be the case even though
echo the program runs fine.
echo " "

# perl waits for 1 second every 1000 lines so it should at least take:
# nlines / 1000 / nProc
# = 90k  / 1000 / 10     = 90 sec

sdfMultiplexer.pl -in .sdf -out t.sdf \
   -cmd "perl -pe '"'$i++ % 1000 == 0 && sleep(1)'"'" -nProc 10 \
>& timetest.err
fgrep -c '$$$$' t.sdf
rm t.sdf

# check to see if less than 40 sec elapsed since we started
if( { (\ps -p $pid -o etime= | egrep -q '^ *00:[0123][0-9]') } ) then
   echo sdfMultiplexer -nProc 10 is significantly faster than single CPU
else
   echo should have taken significantly less than 40sec but was:
       \ps -p $pid -o etime
endif

