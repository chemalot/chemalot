#!/usr/bin/env perl
use warnings;

$script=$0;
$script =~ /^\// || ( $script="$ENV{PWD}/$script" );
our($installDir) = $script =~ /(.*)\/[^\/]+$/;

use File::Copy;
use Getopt::Long;
use POSIX ":sys_wait_h";

my($omegaOpt ) = " -ewindow 20 -maxconfs 500 ";
my($szybkiOpt) = "-MMFF94S -max_iter 2000 -sheffield -grad_conv 1e-8";

$use = <<USE;
sdfMMConfAnalysis.pl [-ref alignRef.sdf] [-szybkiOpt s] [-omegaOpt s] [-optH]
                      -in in.sdf -out out.sdf [-bad bad.sdf] [-sampleOtherMin]
                      [-outputBest maxRMSD]
  -optH .............. optimize hydrogen's before analysis
  -ref  .............. align conformers to given 3d sss fragment
  -szybkiOpt ......... additional option for szybki. default "$szybkiOpt"
                       local optimization will allways use -conj
  -omegaOpt .......... additional option for omega. def "$omegaOpt"
  -sampleOtherMin .... also look for other low energy low rmsd minima
  -harmonic .......... Use harmonic constraints instead of flat buttom
  -failOnMultiInput .. Fail if more than on input record. (workaround for moe bug)
  -bad ............... if given compounds with problems will go to bad.sdf
  -nCPU n ............ Multiplex szibki of omega. (For big inputs lists use
                       multiplex on input)
  -outputBest ........ Output lowest energy pose with rmsd <= maxRMSD
  -debug ............. print commands before executing

If used with sdfMultiplexer use the -groupByAtomCount to make sure the records 
for one input structre are kept together.
USE

my( $com, $in, $ref, $out, $tmpBase, $optH, $debug, $nCPU, $sampleOtherMin,
    $harmonic, $failOnMultiInput, $maxRMSD ) 
  = ("", "", "", "", "", "", "", "", "", "", "", 999);
GetOptions("in=s" =>\$in, 
           "out=s" =>\$out,
           "bad=s" =>\$bad,
           "ref=s" =>\$ref,
           "nCPU=s" =>\$nCPU,
           "optH" => \$optH,
           "outputBest=f" => \$maxRMSD,
           "debug" => \$debug, 
           "harmonic" => \$harmonic, 
           "sampleOtherMin" => \$sampleOtherMin,
           "failOnMultiInput" => \$failOnMultiInput,
           "szybkiOpt=s" => \$szybkiOpt, "omegaOpt=s" => \$omegaOpt) || die $use;
$#ARGV == -1 || die $use;

$in && $out || die $use;
$tmpBase = "/tmp/mmca.$$";
$ref || ($ref = "$tmpBase.sdf");
$maxRMSD = $maxRMSD * 1;

my( @constrVals1 ) = ("5", "5", "5", "5", "5");
my( @constrVals2 ) = ("0.2", "0.6", "1.0", "1.4", "1.8");
my( @constrLbls )  = ("0.2", "0.6", "1.0", "1.4", "1.8");
if( $harmonic )
{  @constrVals1 = ("4", "1", ".25", ".0625", ".016");
   @constrVals2 = ("0", "0", "0", "0", "0");
   @constrLbls  = ("1/4", "1", "4", "16", "64");
}

our($sToKcal) = 0.000238845 * 300; # 300 K
*OUT = *STDOUT;
if( $out !~ /^.sdf$/i )
{  open(OUT, ">$out") || die $!;
}

my( $badOut ) = OUT;
if( $bad )
{   open(BAD, ">$bad") || die $!;
    $badOut = BAD;
}

my($entropy) = "-entropy AN";
$szybkiOpt =~ /solventPB/ && ($entropy = "" ); # not supported
our( $omegaSzybki ) = "szybki $szybkiOpt -in .sdf -out .sdf";

open(IN, "sdfTagTool.csh -in $in -out .sdf -remove 'type|deltaG|deltaE|deltaS|inRMSD|slack'|") || die "$in $!";

my($count) = 0;

while($_ = <IN>)
{  $rec .= $_;
   if( /\$\$\$\$/ )
   {  if( $failOnMultiInput && $count )
      {  die "-failOnMultiInput was given and more than one struct entered\n";
      }
      open(TMP, ">$tmpBase.sdf") || die $!;
      print TMP $rec;
      close(TMP);

      my( $ok ) = &mmConfAnalysis($tmpBase, $ref, $optH, $maxRMSD);

      my( $myOut ) = OUT;
      if( ! $ok )
      {  $myOut = $badOut;
      }

      open( TMP, "$tmpBase.o.sdf" ) || die $!;
      while( $_ = <TMP> ) 
      {  print $myOut $_;;
      }
      close(TMP);

      $count++;
      $rec = "";
   }
}
close(IN);


if( ! $debug )
{  unlink( "$tmpBase.no.sdf", "$tmpBase.hc.1.sdf", "$tmpBase.hc.2.sdf",
           "$tmpBase.hc.3.sdf", "$tmpBase.hc.4.sdf", "$tmpBase.hc.5.sdf",
           "$tmpBase.lo.sdf", "$tmpBase.gmin.sdf", "$tmpBase.sdf", "$tmpBase.other.sdf",
           "$tmpBase.o.sdf", "$tmpBase.hopt.sdf", "$tmpBase.omega.sdf",
           "omega2.log", "omega2.parm", "omega2.rpt", "omega2_status.txt",
           "omega2.pending.ism", 
           "szybki.log", "szybki.param", "szybki.status" );
}


if( $bad && -z $bad ) { unlink( $bad ); }

sub mmConfAnalysis
{  my($tmpBase, $ref, $optH, $maxRMSD) = @_;
   my($inSuffix) = "";
   my($isOK) = 1;
   my($omegaPID) = 0;
   $optH && ($inSuffix = ".hopt");

   # find global minimum
   $com = <<COMS;
     sdfConformerSampler.csh -in .sdf -out .sdf -maxConfs 5 \\
     | $omegaSzybki \\
     | sdfAlign.pl -in .sdf -out .sdf -ref $ref -rmsdTag inRMSD -method sss -mirror
COMS
   chomp($com);

   if( $nCPU > 1 )
   {  $com =~ s/\\?\n//g;
      $com = "sdfMultiplexer.pl -in .sdf -out .sdf -nProc $nCPU -cmd '$com'";
   }
   # find global minimum:
   # omega
   # sdfConformerSampler     sample for OH and NH rotors
   # szybki
   # TODO: secondary sort by RMSD to brak ties in energy
   $com = <<COMS;
   babel -quiet -in $tmpBase.sdf -out .sdf -hydrogens add \\
     | szybki -grad_conv 1 -in .sdf -out .sdf \\
     | omega2 $omegaOpt -in .sdf -out .sdf \\
     | $com \\
     | sdfSorter.csh -in .sdf -out .sdf -numeric -sortTag Total_energy \\
     | tee $tmpBase.omega.sdf \\
     | sdfSplicer.csh -in .sdf -out .sdf -count 1 -readAll \\
     | sdfTagTool.csh -in .sdf -out .sdf -remove "VibRot entropy|Translational entropy|Configurational entropy|Total entropy|Total_energy|MMFF VdW|MMFF Coulomb|MMFF Bond|MMFF Bend|MMFF StretchBend|MMFF Torsion|MMFF Improper Torsion|Ligand MMFF Intramol. Energy" \\
     | szybki $szybkiOpt -in .sdf -out .sdf $entropy \\
     > $tmpBase.gmin.sdf
COMS


   $debug && warn "\n$com\n";
   unless($omegaPID = fork())   # run in background
   {  exec($com);
   }
   if( $nCPU <= 1 ) 
   {  waitpid($omegaPID,0); 
      if( $? != 0 ) { die "Error in executing omega subprocess"; }
   };


   # can not use ! in smarts due to szybki sh conflicts
   $com = <<COMS;
   szybki $szybkiOpt -in $tmpBase.sdf -out .sdf -optGeometry none \\
     |sdfAlign.pl -method sss -in .sdf -ref $ref -out .sdf -rmsdTag inRMSD >$tmpBase.no.sdf

   szybki $szybkiOpt -in $tmpBase.sdf -out .sdf -optGeometry Honly \\
     |sdfConformerSampler.csh -in .sdf -out .sdf -maxConfs 20 \\
     |szybki $szybkiOpt -in .sdf -out .sdf -optGeometry Honly \\
     |sdfSorter.csh -in .sdf -out .sdf -numeric -sortTag Total_energy \\
     |sdfSplicer.csh -in .sdf -out .sdf -count 1 \\
     |sdfAlign.pl -method sss -in .sdf -ref $ref -out .sdf -rmsdTag inRMSD >$tmpBase.hopt.sdf

   szybki  $szybkiOpt -in $tmpBase$inSuffix.sdf -out .sdf $entropy \\
     |sdfAlign.pl -method sss -in .sdf -ref $ref -out .sdf -rmsdTag inRMSD >$tmpBase.lo.sdf 

   szybki  $szybkiOpt -in $tmpBase$inSuffix.sdf -out .sdf -harm_constr1 $constrVals1[0] -harm_constr2 $constrVals2[0] \\
     |sdfAlign.pl -method sss -in .sdf -ref $ref -out .sdf -rmsdTag inRMSD >$tmpBase.hc.1.sdf 

   szybki  $szybkiOpt -in $tmpBase$inSuffix.sdf -out .sdf -harm_constr1 $constrVals1[1] -harm_constr2 $constrVals2[1] \\
     |sdfAlign.pl -method sss -in .sdf -ref $ref -out .sdf -rmsdTag inRMSD >$tmpBase.hc.2.sdf

   szybki  $szybkiOpt -in $tmpBase$inSuffix.sdf -out .sdf -harm_constr1 $constrVals1[2] -harm_constr2 $constrVals2[2] \\
     |sdfAlign.pl -method sss -in .sdf -ref $ref -out .sdf -rmsdTag inRMSD >$tmpBase.hc.3.sdf 

   szybki  $szybkiOpt -in $tmpBase$inSuffix.sdf -out .sdf -harm_constr1 $constrVals1[3] -harm_constr2 $constrVals2[3] \\
     |sdfAlign.pl -method sss -in .sdf -ref $ref -out .sdf -rmsdTag inRMSD >$tmpBase.hc.4.sdf

#   szybki  $szybkiOpt -in $tmpBase$inSuffix.sdf -out .sdf -harm_constr1 $constrVals1[4] -harm_constr2 $constrVals2[4] \\
#     |sdfAlign.pl -method sss -in .sdf -ref $ref -out .sdf -rmsdTag inRMSD >$tmpBase.hc.5.sdf
COMS
   $debug && warn "\n$com\n";
   system($com) == 0 || die "Error in szybki\n";

   if( $nCPU > 1 ) 
   {  waitpid($omegaPID,0); # wait for omega
      if( $? != 0 ) { die "Error in executing omega subprocess"; }
   }

   ($eGMin,$sGMin,$rmsdGMin) =
       split(/\t/,`sdf2Tab.csh -suppressHead -in $tmpBase.gmin.sdf -tags 'Total_energy|Total entropy|inRMSD'`);
   ($eLMin,$sLMin,$rmsdLMin) = 
       split(/\t/,`sdf2Tab.csh -suppressHead -in $tmpBase.lo.sdf -tags 'Total_energy|Total entropy|inRMSD'`);
   $rmsdGMin =~ s/\s//g;
   $rmsdLMin =~ s/\s//g;

   # if local minimu lower in energy than "global" minum.
   # replace global min with local unless -debug or -bad
   if( (!$debug || !$bad )
      && (! -e "$tmpBase.gmin.sdf" || -z "$tmpBase.gmin.sdf" || $eLMin + .5 < $eGMin) )
   {  warn("\nProblem with golbal minimum: eLMin=$eLMin, eGMin=$eGMin\n\n");
      copy("$tmpBase.lo.sdf", "$tmpBase.gmin.sdf");
      $eGMin = $eLMin;
      $sGMin = $sLMin;
      $rmsdGMin = $rmsdLMin;

      $isOK = 0;
   }

   # get set of minima which are better than any of the constraint minima
   # either by energy or by rmsd
   my( $grvyCom ) = <<COMS;
      import groovy.transform.Field;
      \@Field double minRMSD = $rmsdGMin;
      \@Field double minE = $eGMin;

       BigDecimal rmsd = f(\$inRMSD);
       BigDecimal e    = f(\$Total_energy);
       if( minRMSD - rmsd >= 0.05 && e - minE + minRMSD - rmsd >= 0.2) 
       /*if( rmsd < minRMSD && e - minE + Math.abs(rmsd - minRMSD) > 0.2) */
       {  minRMSD = rmsd;
          minE = e;
          return true; 
       };
       return false;
COMS
   $grvyCom =~ s/\n//g;

   $com = <<COMS;
   cat $tmpBase.no.sdf $tmpBase.hopt.sdf \\
         $tmpBase.hc.1.sdf $tmpBase.hc.2.sdf \\
         $tmpBase.hc.3.sdf $tmpBase.hc.4.sdf \\
         $tmpBase.lo.sdf $tmpBase.gmin.sdf \\
   > $tmpBase.o.sdf
   # also pass in constraint minimization so that other minima are also
   # either lower in energy or lower in rmsd than any other result
   cat $tmpBase.o.sdf $tmpBase.omega.sdf \\
     | sdfSorter.csh -in .sdf -out .sdf -numeric -sortTag Total_energy -sortTag inRMSD\\
     | sdfGroovy.csh -in .sdf -out .sdf -c '$grvyCom' \\
     | sdfRMSDSphereExclusion.csh -in .sdf -out .sdf -radius .3 -refFile $tmpBase.o.sdf \\
     | sdfSorter.csh -in .sdf -out .sdf -numeric -sortTag inRMSD \\
     | sdfTagTool.csh -in .sdf -out .sdf -counterTag cntr$$ -addCounter \\
                      -format "type=oth Min{cntr$$}" \\
     | sdfTagTool.csh -in .sdf -out .sdf -remove "VibRot entropy|Translational entropy|Configurational entropy|Total entropy|Total_energy|MMFF VdW|MMFF Coulomb|MMFF Bond|MMFF Bend|MMFF StretchBend|MMFF Torsion|MMFF Improper Torsion|Ligand MMFF Intramol. Energy" \\
     | szybki $szybkiOpt -in .sdf -out .sdf $entropy \\
     | sdfTagTool.csh -in .sdf -out .sdf -remove cntr$$ \\
     > $tmpBase.other.sdf
COMS
   if( $sampleOtherMin )
   {  if( $debug )
      {  warn "\n$com\n";
      }
      system($com) == 0 || die "Error in szybki\n";
   }
   
   &appendDeltaE("$tmpBase.no.sdf",    0,   "input", $eGMin);
   &appendDeltaE("$tmpBase.hopt.sdf", .1,   "H opt", $eGMin);
   &appendDeltaE("$tmpBase.hc.1.sdf", $constrLbls[0], "cstr $constrLbls[0]A", $eGMin);
   &appendDeltaE("$tmpBase.hc.2.sdf", $constrLbls[1], "cstr $constrLbls[1]A", $eGMin);
   &appendDeltaE("$tmpBase.hc.3.sdf", $constrLbls[2], "cstr $constrLbls[2]A", $eGMin);
   &appendDeltaE("$tmpBase.hc.4.sdf", $constrLbls[3], "cstr $constrLbls[3]A", $eGMin);
#   &appendDeltaE("$tmpBase.hc.5.sdf", $constrLbls[4], "cstr $constrLbls[4]A", $eGMin);
   &appendDeltaE("$tmpBase.lo.sdf",    5, "loc Min", $eGMin, $sGMin);
   &appendDeltaE("$tmpBase.gmin.sdf", 10, "glb Min", $eGMin, $sGMin);

   $sampleOtherMin && &appendDeltaE("$tmpBase.other.sdf", 9, "", $eGMin, $sGMin);
   
   $com = "";
   if( $maxRMSD < 999 )
   { # output only lowest e pose with rmsd<maxRMSD 
     $com = "|sdfSorter.csh -in .sdf -out .sdf -numeric -sortTag deltaE"
           ."|sdfGroovy.csh -in .sdf -out .sdf -c 'return (f(\$inRMSD)<=$maxRMSD)'"
           ."|sdfSplicer.csh -in .sdf -out .sdf -count 1";
   }
   $com = <<COMS;
     cat $tmpBase.no.sdf $tmpBase.hopt.sdf \\
         $tmpBase.hc.1.sdf $tmpBase.hc.2.sdf \\
         $tmpBase.hc.3.sdf $tmpBase.hc.4.sdf \\
         $tmpBase.lo.sdf $tmpBase.gmin.sdf $tmpBase.other.sdf \\
     |sdfTagTool.csh -in .sdf -out .sdf \\
         -reorder 'type|inRMSD|deltaE|deltaG|deltaS' \\
         -format 'TITLE={type} dE={deltaE} {inRMSD}' \\
     $com >$tmpBase.o.sdf
COMS
   $debug && warn "\n$com\n";
   system($com) == 0 || die "Error in summarizing\n";

   return $isOK;
}



sub appendDeltaE
{  my( $fName, $slack, $type, $eMin, $sMin) = @_;
   my( $sTxt, $eT, $sT ) = ("", 0, 0);

   open( TMPIN, $fName ) || die "$! ($fName)";

   while( $_ = <TMPIN> )
   {  if( /\$\$\$\$/ )
      {  $sTxt .= getDelatE( $eT, $sT, $slack, $type, $eMin, $sMin); 
         $sTxt .= "\$\$\$\$\n";
         ( $eT, $sT ) = (0, 0);
         next;
      }

      $sTxt .= $_;

      if( /<Total_energy>/ )
      {  $_ = <TMPIN>; 
         $sTxt .= $_;
         chomp;
         $eT = $_;
      }
      if( /<Total entropy>/ )
      {  $_ = <TMPIN>;
         $sTxt .= $_;
         chomp;
         $sT = $_;
      }
   }
   close(TMPIN);

   open(TMPIN, ">$fName" ) || die $!;
   print TMPIN $sTxt;
   close(TMPIN);
}



sub getDelatE
{  my( $eT, $sT, $slack, $type, $eMin, $sMin) = @_;
   my( $fTxt ) = "";
 
   $type && ( $fTxt  = sprintf("> <type>\n%s\n\n",$type));
   $fTxt .= sprintf("> <deltaE>\n%0.2f\n\n", ($eT-$eMin));
   if( $sT && $sMin )
   {  $fTxt .= sprintf("> <deltaS>\n%0.2f\n\n", ($sT-$sMin));
      $fTxt .= sprintf("> <deltaG>\n%0.2f\n\n", ($eT-$eMin - $sToKcal * ($sT-$sMin)));
   }
   $fTxt .= sprintf("> <slack>\n%s\n\n", $slack);

   return $fTxt;
}
