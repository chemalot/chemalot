#!/usr/bin/env perl
use warnings;
use Cwd;

my($script)=$0;
$script =~ /^\// || ( $script="$ENV{PWD}/$script" );
my($installDir) = $script =~ /(.*)\/[^\/]+$/;
my($guasTemplateDir) = "$installDir/data/gaussian";

use Getopt::Long;
$use = "gTorsionMultiplexer.pl -qTemplate template.g -bondFile bd -in inputSdf [-nCPU n]\n"
      ."                       [-core fn] -nSteps n [-startTorsion deg] -torsionIncrement deg\n"
      ."                       [-workDir d] [-out fn] [-queue q] [-minimize -constraintStrength s]\n" 
      ."                       [-maxConfsPerStep] [-jobId id] -finalizeScript\n"
      ."  -qTemplate . template of gaussian input containing #XYZ# and #FName#\n"
      ."  -bondFile .. file with 3d coordinates of 4 torsion atoms, must overly with input\n"
      ."  -out fn .... an output filename. if '.sdf' this will wait for completion.\n"
      ."  -jobId ..... Unique id for queuing (def: id from sequence)\n"
      ."  -workDir ... directory for tmp and out files (def: jobId)\n"
      ."               workDir may not contain files with .g extentions\n"
      ."  -queue ..... lsf queue used (def. medium)\n"
      ."  -nCPU ...... number of cpus to reserve and use in gaussian per optimization (1).\n"
      ."  -in <arg>    input file oe-supported, should be 3d\n"
      ."  -core fn ... align to this substructure, default autogenerate core\n"
      ."  -nSteps      Number of torsion scan steps (def. 20)\n"
      ."  -startTorsion ... start torsion angle for scan (def: no change)\n"
      ."  -torsionIncrement increment from one step to next (def. 18)\n"
      ."  -minimize .. Minimize at each step with MMFF\n"
      ."  -maxConfsPerStep  Generate this number of conformers around other rotatable bondsin\n"
      ."               each step. If -minimize is given only the lowest energy conf is generated\n"
      ."  -constraintStrength One of strong (90),medium (45), weak(20), none or a floating point\n"
      ."               number specifying the strength of tethered constrains for -doMinimize\n"
      ."  -eMail .......... send an eMail to this adress when completed\n"
      ."  -finalizeScript . if given this script is called after the execution has completed.\n"
      ."  -debug .......... if given sub commands are pritnned to stderr\n"
      ."\n"
      ."\tAvailable default templates (check: $guasTemplateDir):\n";

while( <$guasTemplateDir/*.g>)
{  s/.*\///;
   $use .= "\t\t$_\n";
}


my( $template, $bondFile,$in, $finalizeScript, $out, $core, $workDir, $jobId) = ('', '', '', '', '', '', '', '');
my( $nSteps, $startTorsion, $torsionIncrement, $queue, $eMail, $debug ) = (20, '', 18, '', '', '' );
my( $nCPU, $minimize, $constraintStrength, $maxConfsPerStep ) = ( 1, "", "", "" );
GetOptions("qTemplate=s" =>\$template,
           "bondFile=s" => \$bondFile,
           "in=s" => \$in,
           "core=s" => \$core,
           "nSteps=i" => \$nSteps,
           "startTorsion=f" => \$startTorsion,
           "torsionIncrement=f" => \$torsionIncrement,
           "minimize" => \$minimize,
           "maxConfsPerStep=i" => \$maxConfsPerStep,
           "constraintStrength=s" => \$constraintStrength,
           "out=s" => \$out,
           "workDir=s" => \$workDir,
           "jobId=s" => \$jobId,
           "queue=s" => \$queue,
           "nCPU=i" => \$nCPU,
           "eMail=s" => \$eMail,
           "debug" => \$debug,
           "finalizeScript=s" => \$finalizeScript
  ) || die "$!\n$use";

$in       || die "-in is required!\n$use";
$bondFile || die "-bondFile is required!\n$use";
$template || die "-qTemplate is required!\n$use";
$queue    || die "-queue is required!\n$use";
$nSteps           || die "-nSteps is required!\n$use";
$torsionIncrement || die "-torsionIncrement is required!\n$use";

$minimize && ($minimize = "-minimize");
$constraintStrength && ($constraintStrength = "-constraintStrength $constraintStrength");
$maxConfsPerStep && ($maxConfsPerStep = "-maxConfsPerStep $maxConfsPerStep");

-e $template || ($template = "$guasTemplateDir/$template");
-e $template || die "Template: $template not found\n$use";

if( $jobId )
{  $jobId =~ /\// && die "jobId may not contain '/': $jobId";
} else
{  $jobId = "ts_".`tabExport.pl tools getJobId -noHeader`;
   chomp($jobId);
}
my($jIdLen) = length($jobId);

my($oJobNameCom) = "";
if( $out !~ /^.sdf(.gz)?$/i )
{  $oJobNameCom = $out;
   $oJobNameCom =~ s/.*\///;
   $oJobNameCom =~ s/\.sdf(.gz)?$//i;
   $oJobNameCom = "|jobName=$oJobNameCom";
}

$workDir || ($workDir = $jobId);
-e $workDir || mkdir($workDir,0770);
@_=<$workDir/*.g>;
if( $#_ >= 0 ) { die("working directory contains *.g file\n");}
@_=<$workDir/*.out>;
if( $#_ >= 0 ) { die("working directory contains *.out file\n");}


# if core file is given align by it, if not generate core in temp file
my( $outSdf ) = "$workDir/output.sdf";
my($coreOpt ) = ( '' );
if( $core )
{  $core = "| sdfAlign.pl -in .sdf -out .sdf -method SSS -ref '$core'";
} else
{  $core = "| sdfAlign.pl -in .sdf -out .sdf -method SSS -ref '$workDir/core.$$.sdf'";
   $coreOpt = "-core $workDir/core.$$.sdf";
}

# add unique id to input records
# also generate a tab separated file with the original sdf fields to merge back
# when gaussain results are created in sdf file
$com = "sdfTagTool.csh -in $in  -out .sdf -counterTag torMolNum -addCounter -format 'torMolNum=${jobId}_{torMolNum}' "
      ."> '$workDir/in.sdf' ";
$debug && warn( "$com\n" );
system($com) == 0 || die "Error reading input $!\n";


# genereate gaussian inputs
$startTorsion ne "" && ($startTorsion = "-startTorsion $startTorsion");

$com = "QTorsionProfileGenerator.csh -in '$workDir/in.sdf' $coreOpt -template '$template' "
      ." -bondFile '$bondFile' -workDir '$workDir' -outNameTag 'torMolNum' -nCPU $nCPU"
      ." -nSteps $nSteps $startTorsion -torsionIncrement $torsionIncrement "
      ." $minimize $constraintStrength $maxConfsPerStep";

$debug && warn( "$com\n" );
system($com) == 0 || die "Error creating gaussian input $!\n";


# submit gaussian jobs
# do this in the background for faster response
my($nError) = 0;
my($nMol) = 0;
foreach my $gFile ( <$workDir/*.g>)
{  my($gBase) = $gFile =~ /([^\/]+)\.g$/;
   $com = "cd '$workDir';bsub -q $queue -J '$jobId.$gBase' -oo '$gBase.qlog' "
         ."-n $nCPU -R 'span[hosts=1]' gOpt.py -in '$gBase.g' &";
   system($com) == 0 || die "Error executing gaussian on '$gBase.g'\n";

   # do not run more than 8 in parralel
   if( ++$nMol > 8)
   {  if( wait > -1 ) { $! == 0 || $nError++; }
   }
}
$nError += &waitAll;
if( $nError == $nMol ) { die "Failed submittting molecules!" }
if( $nError > 0 )      { warn "$nError errors while submitting molecules!" }

# wait for all bsub commands to submit jobs to queue
my($qWait) = "";
$out =~ /^.sdf(.gz)?$/i && ($qWait = "-K");

my($cwd ) = cwd();
my $unix_file = $cwd . "/" . $out;
my $unix_template_file = $cwd . "/" . "/torScan.template.vortexgz";
my $vfs_file = $jobId . ".vfs";

if( $eMail )
{  my $vortex_url = "http://research.gene.com/vortexweb/vortex.jsp";
   my $mac_path =  "/Volumes/smdd/" . substr($cwd, 29) . "/" . $vfs_file;
   my $pc_path = "////resfiles/smdd/" . substr($cwd, 29) . "/" .  $vfs_file;
   $pc_path =~ s/\//\\/g;
   my $mac_url= $vortex_url . "?url1=" . $mac_path;
   my $pc_url = $vortex_url . "?url1=" . $pc_path;
   $eMail = <<EMAIL;
  cat <<CAT|mail -s 'TorsionScan $out completed' -a '$vfs_file' -a '$out' '$eMail'
   Your Torsion scan job '$out' has completed.
   The working directory was: $cwd

   HTTP link to open job in Vortex (Mac): $mac_url

   HTTP link to open job in Vortex (PC) : $pc_url
CAT
EMAIL
}

# submit finilizer job with dependency on all the gaussian jobs
$template =~ s/^.*\///;
$umask = sprintf(" %o",umask());

# if tempalte like hf_mp2 then we are looking for the final single point mp2 energy
$g2SdfOpt = "-outOptimized";
$template =~ /.+_mp2/ && ($g2SdfOpt = "-outFinal -require MP2");

$com = <<COM;
  bsub $qWait -q veryshort -w 'ended("$jobId.*")' -oo '$workDir/$jobId.final.%J.qulog' -J '$jobId.final' <<'BC'
#!/bin/csh -f
    umask $umask
    set sdfList=()
    foreach n ('$workDir'/*.out)
       g2sdf.pl $g2SdfOpt "\$n" > "\$n:r.sdf" &
       set sdfList=(\$sdfList:q "\$n:r.sdf")
    end
    wait
    
    cat \$sdfList:q \\
    | sdfGroovy.csh -in .sdf -out .sdf -c '\$<>fName=fName.replaceAll(".*\\\\/([^/]+)","\\\$1");\$>TITLE=fName; if(\$<>ScanVar_1.equals("-0.0")) ScanVar_1="0.0"; if(ScanVar_1.equals("-180.0")) ScanVar_1="180.0";' \\
    | sdfTagTool.csh -in .sdf -out .sdf -add 'MinMethod=$template$oJobNameCom' -copy TITLE=torMolNum \\
                     -transform 'torMolNum/(.{$jIdLen}_\\\\d+)\\..*/\$1/' \\
    | sdfSorter.csh -in .sdf -out .sdf -sortTag torMolNum -numeric -sortTag ScanVar_1 -sortTag Energy\\
    | sdfAggregator.csh -in .sdf -out .sdf -outputmode all -groupby torMolNum \\
                        -function 'minE=min(Energy)' \\
    | sdfGroovy.csh -in .sdf -out .sdf -c '\$>deltaE=String.format("%.3f",(float)(f(\$Energy)-f(\$minE))*627.509)' \\
    | sdfAggregator.csh -in .sdf -out .sdf -groupby torMolNum -groupby ScanVar_1 \\
                        -function 'deltaE=min(deltaE)' \\
    $core \\
    | sdfSdfMerger.csh  -master .sdf -out .sdf -second '$workDir/in.sdf' -masterTag torMolNum -secondTag torMolNum \\
    > '$outSdf'
    rm \$sdfList:q

    sdf2Tab.csh -in '$outSdf' -tags 'Energy|ScanVar_1|deltaE|torMolNum' >'$workDir/output.tab'
    if( "$qWait" == "" ) sdfTagTool.csh -in '$outSdf' -out '$out'
    $finalizeScript
    zip --junk-paths $vfs_file '$unix_file' $unix_template_file
    $eMail
BC
COM
$debug && warn( "$com\n" );
system($com) == 0 || die "Error submitting aggregation of results\n";

$qWait && system("sdfTagTool.csh -in '$outSdf' -out '$out'");


sub waitAll
{  my($nError) = 0;
   while( wait > -1 )
   { $! == 0 || $nError++;
   }

   return $nError;
}
