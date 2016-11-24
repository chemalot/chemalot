#!/usr/bin/env perl
use warnings;
use strict;

use Getopt::Long;
use File::Copy;
use File::Path;

#Example commands are in ~smdi/dev/bin/R subfolder

our($progName) = "sdfModelCreateValidate.pl\n";
our($use) = <<___help;
sdfModelCreateValidate.pl\n
   -in <arg>    Supported file format: .sdf
   -out <arg>   Supported file format: .sdf
   -modelType <RandomForest|SVN|Cubist>
   -modelName <arg>   A directory with this name will be created in the
         current directory. The directory will contain all modle files and
         can be used with sdfSSSModelCreateor.pl to make predictions.
   -identifierField <arg>   The field/column containing the record
         (compound) identifier. If is not given the creator will
         generated the identifiers.
   -responseField <arg>   The field/column containing the
         to-be-predicted values.
   -descriptorFile <arg> file containg newline separared descriptor names.
   -predictionColName <arg> Name of the output column with the
         prediction. Required if no conversion Script
   -testRepeats <n> number of random teststo run (10)
   -randomSeed <n> seed for random test set creation.
   -creatorOptions Additinonal options to bass to sdf<modelType>Creator
   -h   Output the help page

   Models are created with sdf<modelType>Creator.pl ./<modelName> is used
   as modelLocation. The test models and the validation results are stored in
   subDirectories. E.g. if your model-name is RF1 you will have this structure:
      ./RF1        (contains main model
      ./RF1/test/0 (contains the 0 test model and its validation results.

   The following command can be used to concatenate all results:
       cat <modelName>/<modelName>.out.sdf \
          <modelName>/test/*/trainSet.out.sdf \
          <modelName>/test/*/testSet.out.sdf
   The output will contain the predictions and the following two new columns:
      model            full|train|test
      testIteration    0 - (<testRepeats> - 1)
   The output is similar to what is passed to the -out file however,
   the -out file will have the descriptorcolumns removed.

   The commands used to create the model are in <modelName>/creatModel.csh
   In case of errors look for: <modelName>/*.log and <modelName>/test/*/*.log
___help


my($testRepeats,$creatorOptions) = (2, "");
my($in,$out,$modelType,$modelName,$identifierField,$responseField) = ();
my($descriptorFile,$predictionColName,$randomSeed,$help) = ();

if( ! GetOptions( 'in=s' => \$in,
                  'out=s' => \$out,
                  'modelType=s' => \$modelType, 
                  'modelName=s' => \$modelName,
                  'identifierField=s' => \$identifierField,
                  'responseField=s' => \$responseField,
                  'descriptorFile=s' => \$descriptorFile, 
                  'predictionColName=s' => \$predictionColName,
                  'testRepeats=i' => \$testRepeats,
                  'randomSeed=i' => \$randomSeed,
                  'creatorOptions=s' => \$creatorOptions,
                  'h' => \$help ))
{  die("\n$use");
}

if( $help )
{  warn( "$use");
   exit(0);
}
warn( "$progName \n\n" );

$in                || &exitWithHelp("-in is required");
$out               || &exitWithHelp("-out is required");
$modelType         || &exitWithHelp("-modelType is required");
$modelName         || &exitWithHelp("-modelName is required");
$identifierField   || &exitWithHelp("-identifierField is required");
$responseField     || &exitWithHelp("-responseField is required");
$descriptorFile    || &exitWithHelp("-descriptorFile is required");
$predictionColName || &exitWithHelp("-predictionColName is required");
-e $descriptorFile || &exitWithHelp("$descriptorFile is not found");

-e $modelName && die("Directory $modelName already exists!\n");
mkdir($modelName) || die("Error creating directory $modelName: $!\n");
if($in =~ /^\..+/ )
{  open(OUT, ">$modelName/inputData.sdf") || die "$!";
   while($_=<>) { print(OUT $_); }
   close(OUT);
}else
{  copy($in, "$modelName/inputData.sdf") || die "$!";
}
$in = "$modelName/inputData.sdf";

my($descriptors, $descriptorCount ) = ("", 0);
open(IN, $descriptorFile) || die "$!";
while($_=<IN>)
{  chomp;
   $descriptors .= "$_|";
   $descriptorCount++;
}
close(IN);
$descriptorCount > 0 || die "Too few descriptors: $descriptorCount\n";
chop($descriptors);

if( $creatorOptions !~ /logImportance/ && $modelType eq "RandomForest" )
{  $creatorOptions .= " -logImportance";
}

# we use interactive bsub so that if we kill the script everithing is killed
my($com) = <<'COM';
bsub -I -q veryshort -n 2 -R "rusage[mem=2] span[hosts=1]" \
     -J '#modelName#' >& '#modelName#/bsub.log' <<'bsub_coms' &
#!/bin/tcsh -f

cat '#inputFile#' \
| sdfTagTool.csh -in .sdf -out .sdf -add 'model=full' \
| sdfR#modelType#Creator.pl -in .sdf -out .sdf #creatorOptions# \
    -modelLocation './#modelName#' -modelName '#modelName#' \
    -responseField '#responseField#' -predictionColName '#predCol#' \
    -identifierField '#identifierField#' \
    -descriptorFields '#descriptors#' \
> '#modelName#/#modelName#.out.sdf'
'bsub_coms'
COM

for( my $i = 0; $i<$testRepeats; $i++)
{  mkpath("$modelName/test/$i") || die "Error creating '$modelName/test/$i': $!\n";

   my($comTest) .= <<'___COM';
   bsub -I -q veryshort -n 2 -R "rusage[mem=2] span[hosts=1]" \
     -J '#modelName#.test.#i#' >& '#modelName#/test/#i#/bsub.log' <<'___bsub_coms' &
   #!/bin/tcsh -f

   cat '#inputFile#' \
   | sdfSplicer.csh -in .sdf -out .sdf -rndSeed #randomSeed# -rndFract 0.75 \
       -skipped '#modelName#/test/#i#/testSetIn.sdf' \
   | sdfTagTool.csh -in .sdf -out .sdf -add 'testIteraton=#i#|model=train' \
   | sdfR#modelType#Creator.pl -in .sdf -out .sdf #creatorOptions# \
       -modelLocation './#modelName#/test/#i#' -modelName '#modelName#' \
       -responseField '#responseField#' -predictionColName '#predCol#' \
       -identifierField '#identifierField#' \
       -descriptorFields '#descriptors#' \
   > '#modelName#/test/#i#/trainSet.out.sdf'

   cat '#modelName#/test/#i#/testSetIn.sdf' \
   | sdfTagTool.csh -in .sdf -out .sdf -add 'testIteraton=#i#|model=test' \
   | sdfRModelPredictor.pl -in .sdf -out .sdf \
       -modelLocation './#modelName#/test/#i#' -modelName '#modelName#' \
   > '#modelName#/test/#i#/testSet.out.sdf'
'___bsub_coms'
___COM

   $comTest =~ s/^   ( *)/$1/gm;
   $comTest =~ s/#i#/$i/g;
   my($rSeed) = $randomSeed * $i;
   $randomSeed && ($comTest =~ s/#randomSeed#/$rSeed/g);

   $com = "$com\n\n#####$i $i $i $i\n$comTest";
}

$com =~ s/#inputFile#/$in/g;
$com =~ s/#modelType#/$modelType/g;
$com =~ s/#modelName#/$modelName/g;
$com =~ s/#creatorOptions#/$creatorOptions/g;
$com =~ s/#responseField#/$responseField/g;
$com =~ s/#identifierField#/$identifierField/g;
$com =~ s/#descriptors#/$descriptors/g;
$com =~ s/#predCol#/$predictionColName/g;

open(OUT, ">$modelName/creatModel.csh") || die "$!";
print(OUT "#!/bin/tcsh -f\n\n");
print(OUT "#AutoCreated with $progName\n\n");
print(OUT "$com\n\n");
print(OUT "# wait for all bsubs to complete:\n");
print(OUT "wait\n");
close(OUT);

# create model
system("tcsh -fc \"source '$modelName/creatModel.csh'\"" ) == 0
|| die "Error creating models: $!";


# concatenate all and remove descriptors
$com  = <<COM;
cat $modelName/$modelName.out.sdf \\
    $modelName/test/*/trainSet.out.sdf \\
    $modelName/test/*/testSet.out.sdf \\
| sdfTagTool.csh -in .sdf -out '$out' -remove '$descriptors'
COM
system($com);



sub exitWithHelp
{  #Function variable declaration and assignment to the local variables.
   my( $info ) = @_;
   die  "$info\n\n"
       ."$use";
}

