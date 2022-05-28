#!/usr/bin/env perl

=head1 NAME

  classifier_cv_cmp.pl

=head1 DESCRIPTION

  Using a final database, perform basic 10-fold cross validation and scoring
  using both PECAN & the RDP Classifier.

  1. Split the database into 10 pieces (taxa > 10 sequences only)
  For each classifier:
    - Train the classifier on parts 1-9. 
    - Classify & Time part 10. 
    - Score the classification results using score_classification.pl

  NOTE: This process considers the entire lineage of the correct taxonomy as follows:
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

  classifier_cv_cmp.pl -f <ref-fasta> -a <ref-tax> -l <ref-lineage> -o <out-dir>

=head1 OPTIONS

  "output-dir|o=s"      => \my $outDir,
  "input-fasta|f=s"     => \my $faFile,
  "input-tax|a=s"       => \my $txFile,
  "input-lineage|l=s"   => \my $spLgFile,
  "num-folds|n=i"       => \$nFolds,
  "offset-coef|c=f"     => \$clusterTbl,
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

=item B<--offset-coef, -c2>
  offset = log10(offsetCoef) is the amount of log10 posterior probability below the
  min(pp) that the error threshold is set.  
  For classification. Default value of offset coef: 0.3.

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

  classifier_cv_cmp.pl classifier_cv_cmp.pl -f /usr/local/bin/PECAN/reference_sequences/V4_trimmed_noEuks_nr_Complete.fa -o v4_cv -a /usr/local/bin/PECAN/post_merge_reference_files/V4_spp_new.tx -l /usr/local/bin/PECAN/post_merge_reference_files/V4_spp_new.lineage


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
my $txSizeThld   = 10;
my $cvSpSizeThld = 10;

GetOptions(
  "output-dir|o=s"      => \my $outDir,
  "input-fasta|f=s"     => \my $faFile,
  "startin-tax|a=s"     => \my $txFile,
  "start-lineage|l=s"   => \my $spLgFile,
  "num-folds|n=i"       => \$nFolds,
  "cluster-Tbl|c=s"     => \my $clusterTbl,
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

###########################################################
## Performing stratified (per species) partition of seq IDs
###########################################################

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

########################
## cross-validation loop
########################

my @pCorrectClKnownSpp;
my %clScore;
my @a;

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
  
  ## The training lineage file has to be made otherwise empty .fa files are made
  ## which will cause est_err_thld to fail
  print "\r[$i] Creating training set lineage file                                  ";
  my $trSpLineageFile = "$cvDir/train.lineage";
  $cmd = "select_fullTx.pl -t $trTxFile -f $spLgFile -o $trSpLineageFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  print "\r[$i] Creating seqID to lineage file                                  ";
  my $trainSeqIDLineageFile = "$cvDir/trainsppSeqID.lineage";
  $cmd = "get_lineage.pl -t $trTxFile -l $spLgFile -o $trainSeqIDLineageFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  ## the only difference is that models will be using fewer sequences for modeling big species

  print "\r[$i] Creating seqID to lineage file                                  ";
  my $testSeqIDLineageFile = "$cvDir/testsppSeqID.lineage";
  $cmd = "get_lineage.pl -t $testTxFile -l $spLgFile -o $testSeqIDLineageFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;


  my $mcDir = "$cvDir/mcDir";

#### Re-format PECAN dataset input files for RDP compatibility ##### 
#### Use sppSeqID files -> reformat -> produce RDP taxonomy file -> addfulllineage w/fasta

  print "\n[$i] Producing RDP-compatible files                                   ";
  my $rdpFmtFile = RDPrfmt($trainSeqIDLineageFile, $cvDir);
  my $rdpLin = "$cvDir/ready4train_taxonomy.txt";
  my $rdpFa = "$cvDir/ready4train_seqs.fa";

  print "\n[$i] Producing RDP-compatible fasta file                            ";
  $cmd = "addFullLineage.py $rdpFmtFile $trFaFile > $rdpFa";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  print "\n[$i] Producing RDP-compatible taxonomy file                            ";
  $cmd = "rdp_training_scripts/lineage2taxTrain.py $rdpFmtFile > $rdpLin";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  print "\n[$i] Training RDP classifier                                   ";
  #### Training RDP classifier ##### 
  $cmd = "java -Xmx4g -jar /usr/local/bin/RDPTools/classifier.jar train -s $rdpFa -t $rdpLin -o $cvDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  #### Training RDP classifier ##### 
  $cmd = "cp /usr/local/bin/classifier/samplefiles/rRNAClassifier.properties $cvDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  print "\n[$i] Classifying test set w/RDP classifier                                   ";
  #### Testing RDP classifier ##### 
  my $rdpOut = "$cvDir/testRDPout.txt";
  $cmd = "time java -Xmx4g -jar /usr/local/bin/RDPTools/classifier.jar classify -t $cvDir/rRNAClassifier.properties -o $rdpOut -f allrank $testFaFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  printRDPout($rdpOut, $cvDir);
  my $rdpclass = $cvDir . "/rdp_class_rfmt.txt";

  #### Training SpeciateIT ##### 
  print "\n[$i] Building model tree and creating taxon's reference fasta files      ";
  $cmd = "rm -rf $mcDir; buildModelTree $quietStr -l $trSpLineageFile -i $trFaFile -t $trTxFile -o $mcDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  print "\n[$i] Building MC models                                              ";
  $cmd = "buildMC -t $mcDir/spp_paths.txt -k 8 -d $mcDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  
  if ($ppEmbedding)
  {
    print "\r[$i] Generating log pp tables for internal nodes of the model tree                                              ";
    $cmd = "pp_embedding -d $mcDir";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  }

  print "\n[$i] Estimating error thresholds                                      ";
  $cmd = "est_error_thlds -d $mcDir -c 0.3";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  
  print "\n[$i] Running classify on test fasta file                               ";
  $cmd = "time classify $ppEmbeddingStr -d $mcDir -i $testFaFile -o $cvDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  print "\n[$i] Building MC models                                              ";
  $cmd = "buildMC -t $mcDir/spp_paths.txt -k 8 -d $mcDir -p 0";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  
  if ($ppEmbedding)
  {
    print "\r[$i] Generating log pp tables for internal nodes of the model tree                                              ";
    $cmd = "pp_embedding -d $mcDir";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  }

  print "\n[$i] Estimating error thresholds                                      ";
  $cmd = "est_error_thlds -d $mcDir -c 0.3";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  
  my %trainTx = readTbl($trTxFile);
  my %testTx = readTbl($testTxFile);
  my $clTxFile = "$cvDir/MC_order7_results.txt";
  my %clTx = readTbl($clTxFile);
  my %testLineage = read2colTbl($testSeqIDLineageFile);

  ### PERFORM CLASSIFICATION SCORING FOR TEST SEQUENCES
  my $clIDs;
  my $clusterTbl = $cvDir . "/merge_out/spp_new_clusterTbl.txt";
  my $clScoreFile;
  my $clcmpFile;
  my %clScore;
  my $clstcp;
  my $rdp_score = $cvDir . "/rdp_scores";
  my $si_score = $cvDir . "/si_scores";

  if ($clusterTbl)
  {
    $cmd = "score_classification.pl -p $rdpclass -l $testSeqIDLineageFile -c $clusterTbl -o $rdp_score";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

    $cmd = "score_classification.pl -p $clTxFile -l $testSeqIDLineageFile -c $clusterTbl -o $si_score";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  }
  else
  {
    $cmd = "score_classification.pl -p $rdpclass -l $testSeqIDLineageFile -o $rdp_score";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

    $cmd = "score_classification.pl -p $clTxFile -l $testSeqIDLineageFile -o $si_score";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  }
        
  my $si_clScoreFile = $si_score . "/clScore.txt";
  my $si_clcmpFile = $si_score . "/cl_cmp.txt";
  my $si_Cl = $cvReportsDir ."/si_clScore_[$i].txt";
  my $si_Cmp = $cvReportsDir ."/si_cl_cmp_[$i].txt";
  my $rdp_Out = $cvReportsDir ."/rawRDPout_[$i].txt";
  print "\r[$i] Copying SI test scores to $cvReportsDir                               ";
  $cmd = "cp $si_clScoreFile $si_Cl; cp $si_clcmpFile $si_Cmp; cp $rdpOut $rdp_Out";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  my $rdp_clScoreFile = $rdp_score . "/clScore.txt";
  my $rdp_clcmpFile = $rdp_score . "/cl_cmp.txt";
  my $rdp_Cl = $cvReportsDir ."/rdp_clScore_[$i].txt";
  my $rdp_Cmp = $cvReportsDir ."/rdp_cl_cmp_[$i].txt";
  print "\r[$i] Copying RDP test scores to $cvReportsDir                               ";
  $cmd = "cp $rdp_clScoreFile $rdp_Cl; cp $rdp_clcmpFile $rdp_Cmp";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
}

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

sub RDPrfmt
{
  my $file = shift;
  my $dir = shift;
  my $out = $dir . "/rdp_formatted_lineage.txt";
  open IN, "<$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  open OUT, ">$out" or die "Cannot open $file for reading: $OS_ERROR\n";
  print OUT "SeqID\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n";
  while (<IN>)
  {
    chomp;
    my ($seqID, $g, $f, $o, $c, $p, $d, $sp) = split /\t/,$_;
    $sp =~ s/\r//g;
    print OUT "$seqID\t$d\t$p\t$c\t$o\t$f\t$g\t$sp\n";
  }
  close IN;
  close OUT;
  return $out;
}

sub printRDPout
{
  my $file = shift;
  my $dir = shift;
  my $out = $dir . "/rdp_class_rfmt.txt";
  open IN, "<$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  open OUT, ">$out" or die "Cannot open $file for reading: $OS_ERROR\n";
  while (<IN>)
  {
    chomp;
    my ($seqID, $burn, $r, $rrank, $rc, $d, $drank, $dc, $p, $prank, $pc, $c, $crank, $cc, $o, $orank, $oc, $f, $frank, $fc, $g, $grank, $gc, $sp, $srank, $sc) = split /\t/,$_;
    $sp =~ s/\r//g;
    print OUT "$seqID\t$sp\n";
  }
  close IN;
  close OUT;
  return $out;
}
exit 0;
