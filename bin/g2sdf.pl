#!/usr/bin/env perl
use warnings;
use Getopt::Long;
$use = "g2sdf.pl [-outAll|-outAllJobs|-outOptimized|-outFinal|-outE E] [-calcDeltaE] gau.out... >out.sdf\n"
      ."  Default is -outFinal\n"
      ."  -outE Search and output only record with the exact string-match energy\n"
      ."  -outFinal Ouput the final gemoetry (and energy of each file\n"
      ."  -outOptimized output the geometry of all converged optimizations\n"
      ."  -outAll output each geometry found in the gaussain output file. Note:\n"
      ."      Optimizations using the #T key will not produce geometries for intermediate steps.\n"
      ."  -outAllJobs output final geometry for each job in file\n"
      ."  -require MP2|HF|DFT will only output if the energy is of given type\n"
      ."  -calcDeltaE find the minimum energy and also report relative energies for all entries [kcal/mol]\n"
      ."              this does NOT distinguish entries by smiles!\n"
      ."\n";
my( $outAll, $outAllJobs, $outOptimized, $outFinal,$outE,$calcDeltaE ) = (0, 0, 0, 0, '', '');
our( $eTypeRequired ) = ".*";
GetOptions("outAll" =>\$outAll, "outAllJobs"=> \$outAllJobs, "outOptimized" => \$outOptimized,
           "outFinal" =>\$outFinal,
           "outE=s" =>\$outE,
           "calcDeltaE"=> \$calcDeltaE,
           "require=s" =>\$eTypeRequired ) || die $use;

if( $outAll + $outAllJobs + $outOptimized + $outFinal != 1 && ! $outE )
{  $outFinal = 1;
}

#TODO compute and output SCS MP2 = E(SCS-MP2) = 6/5*(abab) + 1/3*(aaaa + bbbb)
# In Gaussian (or any other program), you need to find the MP2 output 
# and look for the alpha-alpha(aaaa), beta-beta(bbbb) and alpha-beta(abab) energies.

my($lastMol, $lastJob);
my( $lastE ) = "0";
my( $lastEType ) = "";
#my( $stdOriPattern ) = "(Standard|Z-Matrix) orientation:";
my( $stdOriPattern ) = "(Standard|Z-Matrix|Input) orientation:";
my( $inOPattern ) = "Input orientation:";
my( $oStartPattern ) = "";
my( $hfPattern ) = "SCF Done: *E\\(([^)]+)\\) = *([^ ]+) *A\\.U\\.";
my( $mp2Pattern ) = " EUMP2 = *([^ ]+)";
my( $jobPattern ) = "^ #([PTN] | ).*";
my( $initialParamPattern ) = " !    Initial Parameters    !";
my( $optimizedPattern ) = " !   Optimized Parameters   !";
our( $scanFindPattern ) = "! ([^ ]+) .* (Scan|Frozen|frozen)[, ].*!";
our( $scanValuePattern ) = "! [^ ]+ +.\\(([^ ]+)\\) +([^ ]+) +(-DE|Frozen|frozen|Scan).*!";
our( @scanCoords ) = ();
our( $scanValues ) = {};
our( $scanAtoms ) = {};
our( $scanKeyMap ) = {};
our( $maxScanVarNum ) = 1;
my( $celastLine ) = "";
my( $fName ) = "";
my( $lastFName ) = "###";
my( $oBabel );
if(defined $ENV{"OBABLE"}) 
{  $oBabel = $ENV{"OBABLE"}
}else
{  $oBabel = "obabel";
}

if( $calcDeltaE )
{  open(MYOUT, ">/tmp/g2sdf.$$.sdf" ) || die("Could not open tmp file");
   select(MYOUT);
}


while(<>)
{  $fName = $ARGV;

   if( /$initialParamPattern/ )
   {  # look for scan coodiante names
      @scanCoords = ();
      $scanKeyMap = {};
      $scanValues = {};
      $maxScanVarNum = 1;
      &processOptimizedParam();
   }

   # input file changed
   if( $lastFName ne $fName && $lastMol )
   {  $lastMol = &addEToMol($lastMol,$lastE,$lastEType,$lastJob,$lastFName);
      $outOptimized || print($lastMol);
      $lastMol = "";
      $lastE = "0";
      $lastEType = "";
      warn "g2sdf:Starting on $fName\n";
   }
   $lastFName = $fName;

   # new gaussian job
   if( /$jobPattern/ && $lastLine =~ /^ ---+$/ )
   {  if(/^ #[^T].*NoSym/i)
      {  $oStartPat = $stdOriPattern;
      }else
      {  #$oStartPat = $inOPattern;
         $oStartPat = $stdOriPattern;
      }

      if( $outAllJobs && $lastJob && $lastMol )
      {  print &addEToMol($lastMol,$lastE,$lastEType,$lastJob,$fName);
      }
      $lastJob = $_;
      chomp($lastJob);
      if( /opt[(=]/ && ! /opt=\S+restart/i )
      {  # new optimization job
         # and not a restart reset scaned variable recognition
         foreach my $k (keys %scanValues) { delete($scanValues{$k}); }
      }
   }

   if( /$hfPattern/ )
   {  $lastE = $2;
      chomp($lastE);
      $lastEType = $1;
      if( $outE && $lastE eq $outE )
      {  print &addEToMol($lastMol,$lastE,$lastEType,$lastJob,$fName);
      }
   }
   if( /$mp2Pattern/ )
   {  $lastE = $1;
      chomp($lastE);
      $lastEType = 'MP2';
      if( $outE && $lastE eq $outE )
      {  print &addEToMol($lastMol,$lastE,$lastEType,$lastJob,$fName);
      }
   }

   if( /$optimizedPattern/ )
   {  &processOptimizedParam();
      if( $outOptimized )
      {  print &addEToMol($lastMol,$lastE,$lastEType,$lastJob,$fName);
      }
   }

   # found coordiantes section "Orientation"
   if( $oStartPat && /$oStartPat/ )
   {  if( $outAll && $lastMol )
      {  $lastMol = &addEToMol($lastMol,$lastE,$lastEType,$lastJob,$fName);
         print $lastMol;
      }
      $lastMol = &processCoor();
   }
   $lastLine = $_;
}

if( $lastMol && ! $outOptimized && ! $outE ) 
{  $lastMol = &addEToMol($lastMol,$lastE,$lastEType,$lastJob,$fName);
   print $lastMol;
}


if( $calcDeltaE )
{  close(MYOUT) || die;
   select(STDOUT);
   my($com) = <<COMS;
   sdfAggregator.csh -in /tmp/g2sdf.$$.sdf  -outputmode all -out .sdf \\
        -function 'minE=min(Energy)' \\
     | sdfGroovy.csh -in .sdf -out .sdf \\
        -c '\$>deltaE=String.format("%.3f",(float)(f(\$Energy)-f(\$minE))*627.509)'
COMS
   #warn($com);
   system($com);
}


# also called on initial parameters to set eg. @scanCoords
sub processOptimizedParam
{  my($scanC, $countBreakLine ) = ('', 0);
   
   while(<>)
   {  /-{80}/ && ++$countBreakLine == 2 && return;

      # look for scan coordinates from torsion scans in input
      if( /$scanFindPattern/ )
      {  push(@scanCoords,$1);
      }
   
      # check if this is an input card for a scan coordinate and
      # extract its value
      if( /^ ! / )
      {  foreach $scanC (@scanCoords)
         {  if( index($_, " ! $scanC ") == 0 )
            {  my($atoms,$val) = /$scanValuePattern/;
               $scanValues{$scanC} = sprintf("%.1f",$val);
               $scanAtoms{$scanC} = &convertToIdx($atoms);
               if( !$scanKeyMap{$scanC} )
               {  $scanKeyMap{$scanC} = $maxScanVarNum++;
               }
            }
         }
      }
   }
}

sub processCoor
{  my($molBlock,$com);

   $molBlock = <>;
   $molBlock .= <>;
   $molBlock .= <>;
   $molBlock .= <>;
   while(<>)
   {  $molBlock .= $_;
      /--------/ && last;
   }

   return $molBlock;
}

sub addEToMol
{  my($molBlock, $lastE, $lastEType, $lastJob, $fName) = @_;

   if( $lastEType !~ /$eTypeRequired/ ) { return "";}

   $fName =~ s/\.[^.]+$//;

   $com = "$oBabel -ig09 -osdf<<OBABLE"
             ."|sdfTagTool.csh -in .sdf -out .sdf -addSmi\n"
         ." Entering Link 1 = \n"
         ." #T\n"
         ."Standard orientation:\n"
         .$molBlock
         ."Standard orientation:\n"
         ."OBABLE";
   #warn "$com\n";
   $molBlock = `$com`;
   $molBlock =~ s/<SMI>/<SMILES>/;
   $molBlock = "$fName$molBlock";

   $lastE =~ s/D/E/;
   $lastE = sprintf("%.6f", $lastE);
   my( $tags ) = "> <EnergyType>\n$lastEType\n\n"
                ."> <Job>\n$lastJob\n\n" 
                ."> <fName>\n$fName\n\n"
                ."> <Energy>\n$lastE\n\n";
   foreach my $key (keys %scanValues) 
   {  $tags .= "> <ScanVar_$scanKeyMap{$key}>\n$scanValues{$key}\n\n";
      $tags .= "> <ScanAtoms_$scanKeyMap{$key}>\n$scanAtoms{$key}\n\n";
   }
   $molBlock =~ s/^\$\$\$\$/$tags\$\$\$\$/m;

   return $molBlock;
}

sub convertToIdx
{  my(@atList) = split(/,/,$_[0]);
   my($at);
   for $at (@atList) { $at--; }
   return join(" ", @atList);
}
