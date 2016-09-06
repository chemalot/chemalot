#!/usr/bin/perl -w
$script=$0;
$script =~ /^\// || ( $script="$ENV{PWD}/$script" );
our($installDir) = $script =~ /(.*)\/[^\/]+$/;

use Getopt::Long;
$use = "cluster.pl -in .sdf -out .sdf -radius 0.2 [-nCpu n] -- [options]\n"
   ." in and out files may be any openeye types\n"
   ." additional options will be passed to both SphereEx and NN\n"
   ." nCpu defines the number of threads used in the NN step default=1\n"
   ."   eg. -AAPathSim DEFAULT\n"
   ." The following tags are created: centroidIdx, clusterIdx, NNSim\n"
   ." The following tags will be removed if they exist: includeIdx, sphereIdx\n";

my($in, $out, $radius, $nCpu) = ('','','', 1);
GetOptions("in=s" => \$in, "out=s" => \$out, "radius=f" => \$radius,
           "nCpu=i" => \$nCpu );
$simOpt = join(" ", @ARGV);

($in && $out && $radius) || die $use;

#after sdfMCSSSphereExclusion.csh
#    cluster has spehreIdx, includeIdx maxSim
#    others hade spehreIdx maxSim
#after sdfMCSSNNFinder.csh
#   added NNSim, NNIdx, NNId
$com = "sdfTagTool.csh -in $in -out .sdf -remove 'includeIdx|sphereIdx' |"
      ."sdfMCSSSphereExclusion.csh -radius $radius $simOpt -printAll "
      ."-in .sdf -out .sdf";

# split by existance of includeIdx tag centroids -> clu.$$.sdf
#   cluster members to mem.$$.sdf
$com.= "|sdfGroovy.csh -in .sdf -out .sdf -falseOut /tmp/mem.$$.sdf "
      ."-c 'return \$includeIdx.length() > 0'"
      ."|sdfTagTool.csh -in .sdf -out /tmp/clu.$$.sdf -remove 'maxSim'"
      ."  -rename 'includeIdx=centroidIdx|sphereIdx=clusterIdx' -add NNSim=1\n";

# follow up with nearest neighbor search to reassign members to nearest centroid
$com.= "sdfMCSSNNFinder.csh $simOpt -ref /tmp/clu.$$.sdf -in /tmp/mem.$$.sdf "
      ."-nCpu $nCpu -out .sdf -idTag clusterIdx "
      ."|sdfTagTool.csh -in .sdf -out .sdf -remove 'maxSim|NNIdx|sphereIdx' "
      ." -rename 'NNId=clusterIdx' "
      ."|cat /tmp/clu.$$.sdf - ";

# sort by cluster, keeping cetroid first
$com.= "|sdfSorter.csh -in .sdf -numeric -sortTag clusterIdx -desc -sortTag centroidIdx "
      ." -desc -sortTag NNSim -out $out";
#warn($com);
system($com);
unlink("/tmp/clu.$$.sdf", "/tmp/mem.$$.sdf" );
warn("sdfCluster.pl done\n");
