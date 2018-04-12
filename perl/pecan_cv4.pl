#!/usr/bin/env perl

=head1 NAME

  pecan_cv4.pl

=head1 DESCRIPTION

  Perform cross-validation on all bacteria models in given directory. This process considers 
  the entire lineage of the correct taxonomy as follows:
    - Each species has 7 levels of taxonomy, so correct classification to each of those levels
      equals a point of 1.00/8 = 0.125 x the level reached. 
    - Thus correct classification to each level receives a score as follows: 
        Domain (d_):    1 x 0 = 0
        Phylum (p_):    2 x 0.125 = 0.25
        Class (c_):     3 x 0.125 = 0.325
        Order (o_):     4 x 0.125 = 0.5
        Family (f_):    5 x 0.125 = 0.625
        Genus (g_):     6 x 0.125 = 0.75
        Subgenus (sg_): 7 x 0.125 = 0.875
        Species (s_):   8 x 0.125 = 1

        Where 1 is the best score, and 0.125 is the worst score. 

=head1 SYNOPSIS

  pecan_cv4.pl -o <output dir>

=head1 OPTIONS

"input-dir|i=s"       => \my $mcDir,
  "output-dir|o=s"      => \my $outDir,
  "input-fasta|f=s"     => \my $faFile,
  "input-tax|a=s"       => \my $txFile,
  "input-lineage|l=s"   => \my $spLgFile,
  "num-folds|n=i"       => \$nFolds,
  "offset-coef|c=f"     => \$offsetCoef,
  "offset-coef|c2=f"    => \$offsetCoef2,
  "tx-size-thld|t=i"    => \$txSizeThld,
  "cv-sp-size-thld|s=i" => \$cvSpSizeThld,
  "pp-embedding"        => \my $ppEmbedding,
  "skip-err-thld"       => \my $skipErrThld,
  "verbose|v"           => \my $verbose,
  "debug"               => \my $debug,
  "dry-run"             => \my $dryRun,
  "help|h!"             => \my $help,

=over

=item B<--input-dir, -i>
  Input directory containing all.fa, spp.tx, and spp.lineage files. 

=item B<--output-dir, -o>
  Output directory.

=item B<--sequence-file, -f>
  Fasta file.

=item B<--taxonomy-file, -a>
  Sequence IDs => Species annotation

=item B<--lineage-file, -l>
  Species => Full lineage

=item B<--num-folds, -n>
  Number of folds. That is the number of parts into which the data is split.
  Default value: 10

=item B<--offset-coef, -c>
  offset = log10(offsetCoef) is the amount of log10 posterior probability below the
  min(pp) that the error threshold is set.  
  For model cleaning. Default value of offset coef: 0.9.

=item B<--offset-coef, -c2>
  offset = log10(offsetCoef) is the amount of log10 posterior probability below the
  min(pp) that the error threshold is set.  
  For classification. Default value of offset coef: 0.5.

=item B<--tx-size-thld, -t>
  Size of taxon (number of its ref seq's) below which error threshold is computed
  differently than for the taxons of size at or above the threshold.
  Default value: 10.

=item B<--skip-err-thld>
  Apply --skip-err-thld to the classifier

=item B<--cv-sp-size-thld, -s>
  Cross-validation species size threshold. Only species of that size or higher
  will have their ref seq's partitioned. Default value: 10.

=item B<--pp-embedding>
  Run classifier with the --pp-embedding flag.

=item B<--verbatim, -v>
  Prints content of some output files.

=item B<--debug>
  Prints system commands

=item B<--dry-run>
  Print commands to be executed, but do not execute them.

=item B<-h|--help>
  Print help message and exit successfully.

=back


=head1 EXAMPLE

  pecan_cv4.pl --cv-sp-size-thld 50 --offset-coef 0.9 -n 10 -o all_bacteria_CV_April_25_run1

  pecan_cv4.pl --skip-err-thld --cv-sp-size-thld 50 -o all_bacteria_CV_April_25_run2

=cut

use strict;
use warnings;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use List::MoreUtils qw( part );
use List::Util qw( sum );
use Data::Dumper qw(Dumper);
use Cwd;

## use Math::Permute::Array;

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

my $nFolds       = 10;
my $offsetCoef   = 0.9;
my $txSizeThld   = 10;
my $cvSpSizeThld = 10;
my $offsetCoef2 = 0.50;

GetOptions(
  "input-dir|i=s"       => \my $mcDir,
  "output-dir|o=s"      => \my $outDir,
  "input-fasta|f=s"     => \my $faFile,
  "startin-tax|a=s"     => \my $txFile,
  "start-lineage|l=s"   => \my $spLgFile,
  "final-tax|z=s"       => \my $finalTxFile,
  "final-lin|y=s"       => \my $finalLinFile,
  "num-folds|n=i"       => \$nFolds,
  "offset-coef|c=f"     => \$offsetCoef,
  "offset-coef|c2=f"    => \$offsetCoef2,
  "tx-size-thld|t=i"    => \$txSizeThld,
  "cv-sp-size-thld|s=i" => \$cvSpSizeThld,
  "pp-embedding"        => \my $ppEmbedding,
  "skip-err-thld"       => \my $skipErrThld,
  "verbose|v"           => \my $verbose,
  "debug"               => \my $debug,
  "dry-run"             => \my $dryRun,
  "help|h!"             => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);


if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit 0;
}

if ( !$outDir )
{
  warn "\n\n\tERROR: Missing input group name";
  print "\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
# elsif ( !$nFolds )
# {
#   print "\n\n\t10-fold CV\n\n";
# }


my $debugStr = "";
my $quietStr = "--quiet";
if ($debug)
{
  $debugStr = "--debug";
  $quietStr = "";
}

my $verboseStr = "";
if ($verbose)
{
  $verboseStr = "--verbose";
}

my $skipErrThldStr = "";
if ($skipErrThld)
{
  $skipErrThldStr = "--skip-err-thld";
}

my $ppEmbeddingStr = "";
if ($ppEmbedding)
{
  $ppEmbeddingStr = "--pp-embedding"
}

if ( !$mcDir )
{
  my $mcDir = getcwd;
  print "\n\n\tUsing $mcDir for input files";
  print "\n\n";
}

if ( !$spLgFile )
{
  warn "\n\n\tERROR: $spLgFile does not exist";
  print "\n\n";
  exit 1;
}

elsif ( !$faFile )
{
  warn "\n\n\tERROR: $faFile does not exist";
  print "\n\n";
  exit 1;
}

elsif ( !$txFile )
{
  warn "\n\n\tERROR: $txFile does not exist";
  print "\n";
  exit 1;
}

# if ( -l $treeFile )
# {
#   $treeFile = readlink($treeFile);
# }

####################################################################
##                               MAIN
####################################################################

my $startRun = time();
my $endRun = time();
my $runTime = $endRun - $startRun;
my $timeStr;
my $timeMin = int($runTime / 60);
my $timeSec = $runTime % 60;

print "\r--- Creating cross-validation reports directory";
my $cvReportsDir = "$outDir/cv_reports_dir";
my $cmd = "rm -rf $outDir; mkdir -p $cvReportsDir";
print "\tcmd=$cmd\n" if $dryRun || $debug; # || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;


my $offsetStr = sprintf("o%d", int(100*$offsetCoef));
my $cvReportFile = "$cvReportsDir/accuracy_report_$nFolds" . "FCV_$offsetStr" . ".txt";
if ($skipErrThld)
{
  $cvReportFile = "$cvReportsDir/accuracy_report_$nFolds" . "FCV_minSpSize$cvSpSizeThld" . "_skip_error_thld.txt";
}



##
## Performing stratified (per species) partition of seq IDs
##

my %tx = readTbl($txFile);
## Create hash table of species => seqID
my %spTbl;
for (keys %tx)
{
    push @{$spTbl{$tx{$_}}}, $_;
}


my @base;  # This is array of seq IDs of species with less than
           # $cvSpSizeThld. They will always be in the training set.
my @parts; # Array with nFolds elements, each a ref to an array with seq IDs of
           # one of the nFolds parts
my $nBigSpp = 0; # number of species with at least $cvSpSizeThld elements
my %bigSpp;      # $bigSpp{$sp} = size($sp)

## For each species, consider the number of sequences representing it, 
for my $sp (keys %spTbl)
{
  if ( $sp !~ /_sp/ &&  @{$spTbl{$sp}} >= $cvSpSizeThld )
  {
    my @a = @{$spTbl{$sp}};
    $bigSpp{$sp} = @a;
    $nBigSpp++;
    fisher_yates_shuffle(\@a);

    my $i = 0;
    my @spParts = part { $i++ % $nFolds } @a;

    for my $i (0..($nFolds-1))
    {
      push @{$parts[$i]}, @{$spParts[$i]};
    }
  }
  else #if ($sp !~ /_sp/)
  {
    push @base, @{$spTbl{$sp}};
  }
}

if ($nBigSpp==0)
{
  my $nSpp = keys %spTbl;
  print "\n\nNo. all spp:                    $nSpp\n";
  print     "No. spp with at least $cvSpSizeThld seq's: $nBigSpp\n\n";
  print "\n\tNothing to be done. Exiting !!!\n\n";
  exit 0;
}

if (1)
{
  my $nSpp = keys %spTbl;
  print "\n\nNo. all spp:                    $nSpp\n";
  print     "No. spp with at least $cvSpSizeThld seq's: $nBigSpp\n\n";

  my @spp = sort {$bigSpp{$b} <=> $bigSpp{$a}} keys %bigSpp;
  #printFormatedTbl(\%bigSpp, \@spp);
  my $sum = 0;
  for (keys %bigSpp)
  {
    $sum += $bigSpp{$_};
  }

  print "Total number of seq's:      " . ($sum + @base) . "\n";
  print "No. of seq's in big species: $sum\n";

  print "\nNo. parts: " . @parts . "\n";
  for my $p (@parts)
  {
    ##print "Part size: " . (@{$p} + @base) . "\n";
    print "Part size: " . @{$p} . "\n";
  }
  print "\n";
  # print "Seq's in part 1: @{$parts[0]}\n\n";
  # print "size(base): " . @base . "\n";
  # exit 0;
}

##
## cross-validation loop
##
my @pCorrectClKnownSpp;
my %clScore;
my @a;

my $clRepFile= "$cvReportsDir/CV_summary_report.txt";
open REP, ">$clRepFile" or die "Cannot open $clRepFile for appending: $OS_ERROR";
print REP "CV\tTest nSeqs\tSpecies-Level Accuracy\tOverall Accuracy\n";

foreach my $i (0..($nFolds-1))
{
  print "\r[$i] Creating cross-validation directory";
  
  my $cvDir = "$outDir/cvDir";
  my $cmd = "rm -rf $cvDir; mkdir $cvDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug; # || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  ## Create training set of seqIDs

  my @testIDs = @{$parts[$i]};  
  my @trainIDs;

  foreach my $j (0..($nFolds-1))
  {
    if ( $j != $i )
    {
      push @trainIDs, @{$parts[$j]};
    }
  }

  ##push @testIDs, @base; # this would obscure misclassification error on test sequences as @base would be correctly classified
  push @trainIDs, @base;

  my $nTestIDs = @testIDs;
  my $nTrainIDs = @trainIDs;

  ## Print training seqIDs to seqID file 
  print "\r[$i] Writing training set seqIDs to a file                        ";
  my $trIDsFile = "$cvDir/train.seqIDs";
  print_seqID_to_file($trIDsFile, \@trainIDs);

  ## Print testing seqIDs to seqID file 
  print "\r[$i] Writing test seqIDs to a file                                  ";
  my $testIDsFile = "$cvDir/test.seqIDs";
  print_seqID_to_file($testIDsFile, \@testIDs);

  print "\r[$i] Creating training set fasta file                                 ";
  my $trFaFile = "$cvDir/train.fa";
  $cmd = "select_seqs.pl $quietStr -s $trIDsFile -i $faFile -o $trFaFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  print "\r[$i] Creating test set fasta file                                      ";
  my $testFaFile = "$cvDir/test.fa";
  $cmd = "select_seqs.pl $quietStr -s $testIDsFile -i $faFile -o $testFaFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  print "\r[$i] Creating training set taxon file                                  ";
  my $trTxFile = "$cvDir/train.tx";
  $cmd = "select_tx.pl $quietStr -s $trIDsFile -i $txFile -o $trTxFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  print "\r[$i] Creating test taxon file                                           ";
  my $testTxFile = "$cvDir/test.tx";
  $cmd = "select_tx.pl $quietStr -s $testIDsFile -i $txFile -o $testTxFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  my $fullLin = $cvDir . "/full.lineage";
  if ($finalTxFile)
  {
    update_test_tx(\$testTxFile, \$finalTxFile);
    $cmd = "cat $spLgFile $finalLinFile | sort | uniq > $fullLin";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  }
  
  ## The training lineage file has to be made otherwise empty .fa files are made
  ## which will cause est_err_thld to fail
  print "\r[$i] Creating training set lineage file                                  ";
  my $trSpLineageFile = "$cvDir/train.lineage";
  $cmd = "select_fullTx.pl -t $trTxFile -f $spLgFile -o $trSpLineageFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  
  my $mcDir = "$cvDir/mcDir";

  print "\r[$i] Creating seqID to lineage file                                  ";
  my $trainSeqIDLineageFile = "$cvDir/trainsppSeqID.lineage";
  $cmd = "get_lineage.pl -t $trTxFile -l $spLgFile -o $trainSeqIDLineageFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  ## the only difference is that models will be using fewer sequences for modeling big species

  print "\r[$i] Creating seqID to lineage file                                  ";
  my $testSeqIDLineageFile = "$cvDir/testsppSeqID.lineage";
  $cmd = "get_lineage.pl -t $testTxFile -l $fullLin -o $testSeqIDLineageFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  print "\r[$i] Building model tree and creating taxon's reference fasta files      ";
  $cmd = "rm -rf $mcDir; buildModelTree $quietStr -l $trSpLineageFile -i $trFaFile -t $trTxFile -o $mcDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  print "\r[$i] Building MC models                                              ";
  $cmd = "buildMC -t $mcDir/spp_paths.txt -k 8 -d $mcDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  print "\r[$i] Running classify on training fasta file                               ";
  $cmd = "time classify -d $mcDir -i $trFaFile -o $cvDir --skip-err-thld";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  my $mcclass = $cvDir . "/MC_order7_results.txt";
  print "\r[$i] Scoring classification of training fasta file                               ";
  $cmd = "score_classification.pl -l $trainSeqIDLineageFile -p $mcclass -o $cvDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  my $mcScore = $cvDir . "/clScore.txt";
  my $mcCmp = $cvDir . "/cl_cmp.txt";
  my $beforeCl = $cvReportsDir ."/train_clScore_beforeMerge_[$i].txt";
  my $beforeCmp = $cvReportsDir ."/train_cl_cmp_beforeMerge_[$i].txt";

  print "\r[$i] Copying training before scores to $cvReportsDir                               ";
  $cmd = "cp $mcScore $beforeCl; cp $mcCmp $beforeCmp";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  my $mergeOut = $cvDir . "/merge_out";
  print "\r[$i] Merging models                               ";
  print "\r[$i] Writing to $mergeOut                               ";
  $cmd = "merge_models_2.pl -p $mcclass -l $trSpLineageFile -t $trTxFile -c $mcScore -o $mergeOut";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  my $newTx = $cvDir . "/merge_out/spp_new.tx";
  my $newLin = $cvDir . "/merge_out/spp_new.lineage";

  print "\r[$i] Creating seqID to lineage file                                  ";
  $cmd = "get_lineage.pl -t $newTx -l $newLin -o $trainSeqIDLineageFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  print "\r[$i] Building model tree and creating taxon's reference fasta files      ";
  $cmd = "rm -rf $mcDir; buildModelTree $quietStr -l $newLin -i $trFaFile -t $newTx -o $mcDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  print "\r[$i] Building MC models                                              ";
  $cmd = "buildMC -t $mcDir/spp_paths.txt -k 8 -d $mcDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  print "\r[$i] Running classify on training fasta file                               ";
  $cmd = "time classify -d $mcDir -i $trFaFile -o $cvDir --skip-err-thld";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  $mcclass = $cvDir . "/MC_order7_results.txt";
  print "\r[$i] Scoring classification of training fasta file                               ";
  $cmd = "score_classification.pl -l $trainSeqIDLineageFile -p $mcclass -o $cvDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  $mcScore = $cvDir . "/clScore.txt";
  $mcCmp = $cvDir . "/cl_cmp.txt";
  my $afterCl = $cvReportsDir ."/train_clScore_afterMerge_[$i].txt";
  my $afterCmp = $cvReportsDir ."/train_cl_cmp_afterMerge_[$i].txt";

  print "\r[$i] Copying training before scores to $cvReportsDir                               ";
  $cmd = "cp $mcScore $afterCl; cp $mcCmp $afterCmp";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  
  if ($ppEmbedding)
  {
    print "\r[$i] Generating log pp tables for internal nodes of the model tree                                              ";
    $cmd = "pp_embedding -d $mcDir";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  }

  if (!$skipErrThld)
  {
    print "\r[$i] Estimating error thresholds                                      ";
    $cmd = "est_error_thlds -d $mcDir -c $offsetCoef2";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  }
  
  print "\r[$i] Running classify on test fasta file                               ";
  $cmd = "time classify $ppEmbeddingStr $skipErrThldStr -d $mcDir -i $testFaFile -o $cvDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  if ($ppEmbedding)
  {
    $cmd = "cat $mcDir/classification_paths.txt > $cvReportsDir/classification_paths_[$i].txt";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  }

  print "\r[$i] Creating a table of true and classifier's taxonomy                 ";
  my %trainTx = readTbl($trTxFile);
  my %testTx = readTbl($testTxFile);
  my $clTxFile = "$cvDir/MC_order7_results.txt";
  my %clTx = readTbl($clTxFile);
  my %testLineage = read2colTbl($testSeqIDLineageFile);

  if ($debug)
  {
    print Dumper \%testLineage;
  }

  ## Checking if both tables have the same seqIDs as their keys
  my @testTxIDs = keys %testTx;

  if ( @testTxIDs != @testIDs )
  {
    warn "ERROR: size(testTxIDs) != size(testIDs)";
    print "size(testTxIDs): " . @testTxIDs . "\n";
    print "size(testIDs): " . @testIDs . "\n";
    exit 1;
  }

  
  ## Checking if both the classified and test table have the same seqIDs as their keys
  my @clIDs   = keys %clTx;
  my @c = comm(\@testIDs, \@clIDs);
  if (@c != @testIDs || @c != @clIDs)
  {
    warn "WARNING: seq IDs of the testTx and clTx tables do not match";
    print "Number of elements in testTx: " . @testIDs . "\n";
    print "Number of elements in clTx: " . @clIDs . "\n";
    print "Number of elements common to both tables: " . @c . "\n\n";

    if ( @testIDs > @clIDs )
    {
      my @d = diff(\@testIDs, \@clIDs);
      print "\nElements in testTx that are not clTx:\n";
      for (@d)
      {
	     print "\t$_\t" . $testTx{$_} . "\n";
      }
      print "\n\n";
    }

    if ( @clIDs > @testIDs )
    {
      my @d = diff(\@clIDs, \@testIDs);
      print "\nElements in clTx that are not in testTx:\n";
      for (@d)
      {
	     print "\t$_\t" . $clTx{$_} . "\n";
      }
      print "\n\n";
    }
    exit 1;
  }

  ### PERFORM CLASSIFICATION SCORING FOR TEST SEQUENCEScd 
  my $clIDs;
  my $clusterTbl = $cvDir . "/merge_out/spp_new_clusterTbl.txt";
  my $clScoreFile;
  my $clcmpFile;
  my %clScore;
  my $clstcp;
  if ($clusterTbl)
  {
    my $clstcp = $cvReportsDir ."/train_clusterTbl_[$i].txt";
    print "\r[$i] Copying cluster table to $cvReportsDir                               ";
    $cmd = "cp $clusterTbl $clstcp";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

    $cmd = "/usr/local/bin/PECAN/bin/score_classification.pl -p $clTxFile -l $testSeqIDLineageFile -c $clusterTbl -o $cvDir";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
    
    $clScoreFile = $cvDir . "/clScore.txt";
    $clcmpFile = $cvDir . "/cl_cmp.txt";
    $afterCl = $cvReportsDir ."/test_clScore_[$i].txt";
    $afterCmp = $cvReportsDir ."/test_cl_cmp_[$i].txt";
    print "\r[$i] Copying training before scores to $cvReportsDir                               ";
    $cmd = "cp $clScoreFile $afterCl; cp $clcmpFile $afterCmp";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
    
    %clScore = readTbl($clScoreFile);
  }
  else
  {
    $cmd = "/usr/local/bin/PECAN/bin/score_classification.pl -p $clTxFile -l $testSeqIDLineageFile -o $cvDir";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
    
    $clScoreFile = $cvDir . "/clScore.txt";
    $clcmpFile = $cvDir . "/cl_cmp.txt";
    $afterCl = $cvReportsDir ."/test_clScore_[$i].txt";
    $afterCmp = $cvReportsDir ."/test_cl_cmp_[$i].txt";
    print "\r[$i] Copying training before scores to $cvReportsDir                               ";
    $cmd = "cp $clScoreFile $afterCl; cp $clcmpFile $afterCmp";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
    
    %clScore = readTbl($clScoreFile);
  }
  my $clScore_sum = sum values %clScore;
  my $clScore_perc = 100.0 * $clScore_sum /  scalar keys %clScore;

  my @ones = grep {$clScore{$_} eq "1"} keys %clScore;
  my $sppScore = scalar (@ones);
  my $sppclScore_perc = 100.0 * $sppScore /  scalar keys %clScore;

  if ($skipErrThld) 
  {
    print "\n\nWithout error thresholds, the species-level classification score for cv_[$i] is: " . $sppclScore_perc . "\n";
    print "\n\nWithout error thresholds, the overall classification score for cv_[$i] is: " . $clScore_perc . "\n";
  }
  else 
  {
    print "\n\nUsing error thresholds, the species-level classification score for cv_[$i] is: " . $sppclScore_perc . "\n";
    print "\n\nUsing error thresholds, the overall classification score for cv_[$i] is: " . $clScore_perc . "\n";
  }

  print REP "$i\t".length (keys %clScore)."\t$sppclScore_perc\t$clScore_perc\n";
}
close REP;

$endRun = time();
$runTime = $endRun - $startRun;
if ( $runTime > 60 )
{
  $timeMin = int($runTime / 60);
  $timeSec = sprintf("%02d", $runTime % 60);
  $timeStr = "$timeMin:$timeSec";
}
else
{
  $runTime = sprintf("%02d", $runTime);
  $timeStr = "$timeMin:$runTime";
}

$startRun = time();
print "\n\tCompleted in $timeStr\n\n";
print "\tReport written to $cvReportFile\n\n";

####################################################################
##                               SUBS
####################################################################

# read lineage table
sub readLineageTbl{

  my $file = shift;

  if ( ! -f $file )
  {
    warn "\n\n\tERROR in readLineageTbl(): $file does not exist";
    print "\n\n";
    exit 1;
  }

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  foreach (<IN>)
  {
    chomp;
    my ($id, $t) = split /\s+/,$_;
    $tbl{$id} = $t;
    ## test for '/' characters
    if ($t =~ /\//)
    {
      warn "\n\n\tERROR: Discovered '/' for id: $id\t$t";
      print "\n\n";
      exit 1;
    }
  }
  close IN;

  return %tbl;
}

sub get_seqIDs_from_fa
{
  my $file = shift;

  my $quiet = 1;
  my $startRun = time();
  my $endRun = time();

  open (IN, "<$file") or die "Cannot open $file for reading: $OS_ERROR\n";
  $/ = ">";
  my $junkFirstOne = <IN>;
  my $count = 1;
  my $timeStr = "";
  my @seqIDs;
  while (<IN>)
  {
    if ( !$quiet && ($count % 500 == 0) )
    {
      $endRun = time();
      my $runTime = $endRun - $startRun;
      if ( $runTime > 60 )
      {
	my $timeMin = int($runTime / 60);
	my $timeSec = sprintf("%02d", $runTime % 60);
	$timeStr = "$timeMin:$timeSec";
      }
      else
      {
	my $runTime = sprintf("%02d", $runTime);
	$timeStr = "0:$runTime";
      }
      print "\r$timeStr";
    }

    chomp;
    my ($id,@seqlines) = split /\n/, $_;
    push @seqIDs, $id;
    $count++;
  }
  close IN;
  $/ = "\n";

  return @seqIDs;
}

# common part of two arrays
sub comm{

  my ($a1, $a2) = @_;

  my @c; # common array
  my %count;

  foreach my $e (@{$a1}, @{$a2}){ $count{$e}++ }

  foreach my $e (keys %count)
  {
    push @c, $e if $count{$e} == 2;
  }

  return @c;
}

# read table with one column
sub readArray{

  my ($file, $hasHeader) = @_;
  my @rows;

  if ( ! -f $file )
  {
    warn "\n\n\tERROR in readArray(): $file does not exist";
    print "\n\n";
    exit 1;
  }

  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  if ( defined $hasHeader )
  {
    <IN>;
  }
  foreach (<IN>)
  {
    chomp;
    push @rows, $_;
  }
  close IN;

  return @rows;
}

# fisher_yates_shuffle( \@array ) : generate a random permutation
# of @array in place
sub fisher_yates_shuffle {
    my $array = shift;
    my $i;
    for ($i = @$array; --$i; ) {
        my $j = int rand ($i+1);
        next if $i == $j;
        @$array[$i,$j] = @$array[$j,$i];
    }
}

# read two column table; create a table that assigns
# elements of the first column to the second column
sub readTbl{

  my $file = shift;

  if ( ! -f $file )
  {
    warn "\n\n\tERROR in readTbl(): $file does not exist";
    print "\n\n";
    exit 1;
  }

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  foreach (<IN>)
  {
    chomp;
    my ($id, $t) = split /\s+/,$_;
    $tbl{$id} = $t;
  }
  close IN;

  return %tbl;
}

# difference of two arrays
sub diff{

  my ($a1, $a2) = @_;

  my (%aa1, %aa2);

  foreach my $e (@{$a1}){ $aa1{$e} = 1; }
  foreach my $e (@{$a2}){ $aa2{$e} = 1; }

  my @d; # dfference array

  foreach my $e (keys %aa1, keys %aa2)
  {
    push @d, $e if exists $aa1{$e} && !exists $aa2{$e};
  }

  return @d;
}

# extract unique elements from an array
sub unique{

  my $a = shift;
  my %saw;
  my @out = grep(!$saw{$_}++, @{$a});

  return @out;
}

# print elements of a hash table
sub printTbl{

  my $rTbl = shift;
  map {print "$_\t" . $rTbl->{$_} . "\n"} keys %$rTbl;
}

# print elements of a hash table so that arguments are aligned
sub printFormatedTbl{

  my ($rTbl, $rSub) = @_; # the second argument is a subarray of the keys of the table

  my @args;
  if ($rSub)
  {
    @args = @{$rSub};
  }
  else
  {
    @args = keys %{$rTbl};
  }

  my $maxStrLen = 0;
  map { $maxStrLen = length($_) if( length($_) > $maxStrLen )} @args;

  for (@args)
  {
    my $n = $maxStrLen - length($_);
    my $pad = ": ";
    for (my $i=0; $i<$n; $i++)
    {
      $pad .= " ";
    }
    print "$_$pad" . $rTbl->{$_} . "\n";
  }
  print "\n";
}

sub read2colTbl{

  my $file = shift;

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  foreach (<IN>)
  {
    chomp;
    next if $_ eq "";
    my ($id, $t) = split (/\s+/,$_, 2);
    $tbl{$id} = $t;
  }
  close IN;

  return %tbl;
}

sub print_seqID_to_file
{
  my $IDsFile = shift;
  my $rIDs = shift;

  my @IDs = @{ $rIDs};

  open OUT, ">$IDsFile" or die "Cannot open $IDsFile for writing: $OS_ERROR";
  for my $id (@IDs)
  {
    print OUT "$id\n";
  }
  close OUT;
}

sub update_test_tx
{
  my $rtestTxFile = shift;
  my $rfinalTxFile = shift;

  my $testTxFile = $$rtestTxFile;
  my $finalTxFile = $$rfinalTxFile;

  my %tx = read2colTbl($testTxFile);
  my %newtx = read2colTbl($finalTxFile);

  open OUT, ">$testTxFile" or die "Cannot open $testTxFile for writing: $OS_ERROR";
  foreach my $x (keys %tx)
  {
    my $old = $tx{$x};
    my $new = $newtx{$x};

    if ($new =~ /$old/ || $new =~ /Cluster/ )
    {
      print OUT "$x\t$old\n";
    }
    else
    {
      print OUT "$x\t$new\n";
    }
  }
  close OUT;
}
exit 0;
