#!/bin/csh -f

################################################
# To run a simple test and execute most commands issue the following commands from the chemalot directory:
#
# test/simple.csh
################################################

(echo CCF;echo CCO )\
| sdfFingerprinter.csh -in .smi -out .sdf -format folded512 -fpType 'linear7*4' \
| sdfFingerprinter.csh -in .sdf -out .sdf -format folded512 -fpType 'HashLinear7*4' \
| sdfFingerprinter.csh -in .sdf -out .sdf -format folded512 -fpType 'maccs'\
| sdfFPNNFinder.csh -in .sdf -out .sdf -fpTag 'linear7*4_folded512' \
| sdfTagTool.csh -in .sdf -out .sdf -rename NNSim=74.NNSim \
| sdfMCSSNNFinder.csh -in .sdf -out .sdf -AAPathSim DEFAULT \
| sdfTagTool.csh -in .sdf -out .sdf -rename NNSim=AA.NNSim \
| sdfSorter.csh -in .sdf -out .sdf -desc -numeric -sortTag AA.NNSim \
| sdfCluster.pl -in .sdf -out .sdf -radius 0.3 -- -AAPathSim DEFAULT \
| perl -pe 's/^(  -OEChem-).+/$1/' \
>test/simple.out.sdf

# compare actuall output to expected output
diff test/simple.ref.sdf test/simple.out.sdf

