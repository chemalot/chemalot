#!/usr/bin/env perl
use warnings;
$script=$0;
$script =~ /^\// || ( $script="$ENV{PWD}/$script" );
our($installDir) = $script =~ /(.*)\/[^\/]+$/;

use Getopt::Long;
$use = "cluster.pl -in .sdf -out .sdf -radius 0.2 -fpTag tag -- options\n"
   ." in and out files may be any openeye types\n"
   ." additional options will be passed to both SphereEx and NN\n"
   ."\n";

my($in, $out, $radius, $fpTag) = ('','','');
GetOptions("in=s" => \$in, "out=s" => \$out, "radius=f" => \$radius,
           "fpTag=s" => \$fpTag);
$simOpt = join(" ", @ARGV);

($in && $out && $radius && $fpTag) || die $use;

#after sdfFPSphereExclusion.csh
#    cluster has spehreIdx, includeIdx maxSim
#    others hade spehreIdx maxSim
#after sdfFPNNFinder.csh
#   added NNSim, NNIdx, NNId
$com = "sdfFPSphereExclusion.csh -fpTag $fpTag -radius $radius $simOpt -printAll "
      ."-in $in -out .sdf";

# split by existance of includeIdx tag centroids -> clu.$$.sdf
#   cluster members to mem.$$.sdf
$com.= "|sdfGroovy.csh -in .sdf -out .sdf -falseOut /tmp/mem.$$.sdf "
      ."-c 'return \$includeIdx.length() > 0'"
      ."|sdfTagTool.csh -in .sdf -out /tmp/clu.$$.sdf -remove 'maxSim'"
      ."  -rename 'includeIdx=centroidIdx|sphereIdx=clusterIdx' -add NNSim=1\n";

# follow up with nearest neighbor search to reassign members to nearest centroid
$com.= "sdfFPNNFinder.csh -fpTag $fpTag $simOpt -ref /tmp/clu.$$.sdf -in /tmp/mem.$$.sdf "
      ."-out .sdf -idTag clusterIdx "
      ."|sdfTagTool.csh -in .sdf -out .sdf -remove 'maxSim|NNIdx|sphereIdx' "
      ." -rename 'NNId=clusterIdx' "
      ."|cat /tmp/clu.$$.sdf - ";

# sort by cluster, keeping cetroid first
$com.= "|sdfSorter.csh -in .sdf -numeric -sortTag clusterIdx -desc -sortTag centroidIdx "
      ." -desc -sortTag NNSim -out $out";
#warn($com);
system($com);
unlink("/tmp/clu.$$.sdf", "/tmp/mem.$$.sdf" );

