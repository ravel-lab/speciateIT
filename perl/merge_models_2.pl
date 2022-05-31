#!/usr/bin/env perl

=head1 NAME

  merge_models.pl

=head1 DESCRIPTION

  Using the classification scores produced by score_classification.pl, determine which PECAN
  models can be merged by determining the reference species that are misclassified, which species the
  reference sequences are misclassified to, and merging the models together by re-labeling the
  associated taxonomy and altering the species annotation in the lineage file.

  1. For each sequence, determine if the PECAN and reference classification match.
      -> Use clScore.txt to get seqIDs w/score of 0
  2. For each key of the hash of arrays, call the array and for each item, attach it to a single label.
  3. Grep the taxonomy file hash for values equal to the keys of the hash of arrays, and for each of them, re-associated the
      sequence ID to the new, concanated label.
      -> print hash to new taxonomy file. (spp_new.tx)
  4. Grep the lineage hash file for keys equal to the keys of the hash of arrays, and for each of these, re-associate the
      sequence ID to the new, concatenated label.
      -> print hash to new lineage file (spp_new.lineage)

=head1 SYNOPSIS

  merge_models.pl -i <input dir> -o <output dir>

  **Currently, looks for sppSeqID.lineage and MC_order7_results.txt
    as input from the in directory.

  "lineage-file|l=s"    => \my $LineageFile,
  "taxonomy-file|t=s"   => \my $TxFile,
  "class-scores|c=s"    => \my $clScoreFile,
  "pecan-output|p=s"    => \my $pecanFile,
  "offset-coef|z=s"     => \my $offsetCoef,
  "tx-size-thld|n=s"    => \my $txSizeThld,
  "pp-embedding"        => \my $ppEmbedding,
  "skip-err-thld"       => \my $skipErrThld,
  "quiet"               => \my $quiet,
  "verbose|v"           => \my $verbose,
  "debug"               => \my $debug,
  "dry-run"             => \my $dryRun,
  "help|h!"             => \my $help,
  )

=head1 OPTIONS

=over

=item B<--input-dir, -i>
  Input directory containing spp.tx, and spp.lineage, clScore.txt, and MC_order7_results.txt files.

=item B<--output-dir, -o>
  Output directory.

=item B<--lineage-file, -l>
  If a lineage file other than spp.lineage is to be used, specify here.

=item B<--taxonomy-file, -t>
  If a taxonomy file other than spp.tx is to be used, specify here.

=item B<--classification-score-file, -c>
  If a score file other than clScore.txt is to be used, specify here.

=item B<--pecan-output-file, -p>
  If a pecan classification file other than MC_order7_results.txt is to be used, specify here.

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

  improve_models.pl -i ~/Desktop/PECAN/FFT_NS_2/FFT_NS_2/V3V4_models_bac_only -o ~/Desktop/PECAN/FFT_NS_2/FFT_NS_2

=cut

use strict;
use warnings;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use List::MoreUtils qw( part );
use List::Util qw( sum );
use Data::Dumper qw(Dumper);
use List::MoreUtils qw(uniq);
use Cwd;
use Devel::Size qw(size total_size);

## use Math::Permute::Array;

$OUTPUT_AUTOFLUSH = 1;
####################################################################
##                             OPTIONS
####################################################################
GetOptions(
  "lineage-file|l=s"    => \my $LineageFile,
  "taxonomy-file|t=s"   => \my $TxFile,
  "class-scores|c=s"    => \my $clScoreFile,
  "pecan-output|p=s"    => \my $pecanFile,
  "output-dir|o=s"      => \my $outDir,
  "quiet"               => \my $quiet,
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

if ( !$LineageFile && !$TxFile && !$clScoreFile && !$pecanFile)
{
  print "Please provide lineage (-l), taxonomy (-t), clScore (-c), and PECAN (-p) files\n";
  exit 0;
}
elsif ( !$LineageFile )
{
  print "Please provide lineage (-l) file\n";
  exit 0;
}
elsif ( !$TxFile )
{
  print "Please provide taxonomy (-t) files\n";
  exit 0;
}
elsif ( !$clScoreFile )
{
  $clScoreFile = "clScore.txt";
}
elsif ( !$pecanFile)
{
  $pecanFile = "MC_order7_results.txt";
}

# if (!$outDir)
# {
#    $outDir= "merge_out";
# }

my $debugStr;
if ($debug)
{
  $debugStr = "--debug";
}
else
{
  $debugStr = "";
}
print "Perl $^V\n";

####################################################################
##                               MAIN
####################################################################
my $generaDir = $outDir . "/genera_files";
my $cmd;

print "\tUsing $LineageFile for lineage file\n";
print "\tUsing $TxFile for taxonomy file\n";
print "\tUsing $outDir for writing new files\n";

if (-e $outDir)
{
  $cmd = "rm -rf $outDir; mkdir $outDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
}
else
{
  $cmd = "mkdir $outDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
}

$cmd = "rm -rf $generaDir; mkdir $generaDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

## Read in scores, taxonomy, and PECAN output file
my %clScore = readTbl($clScoreFile);
my %pecan = readTbl($pecanFile);
my @misclass;
my $fixedTxFile = $outDir . "/spp_new.tx";
if (-e $fixedTxFile)
{
  unlink $fixedTxFile;
}
my $fixedLinFile = $outDir . "/spp_new.lineage";
if (-e $fixedLinFile)
{
  unlink $fixedLinFile;
}
my %newTx;
my %genera;
my $genus;
my $spp;

foreach my $seqID (sort keys %clScore) ## For each of the sequences
{
  my $score = $clScore{$seqID};
  if ($score == 0)
  {
    push @misclass, $seqID;
    my $taxon = $pecan{$seqID};
    my ($genus, $spp) = split(/_/, $taxon);
    push @{ $genera{$genus} }, $seqID; ## %genera: genus => array of misclassified seqIDs
  }
}


## Launch merging.pl for each genus of misclassified seq's
foreach $genus (sort keys %genera) ## %genera: genus => array of seqIDs
{
  my $misclassFile = "$generaDir/". $genus . "_misclassified.txt";
  open OUT, ">$misclassFile" or die "Cannot open $misclassFile for writing: $OS_ERROR";
  my $rseqs = $genera{$genus};
  my @seqs = @$rseqs;
  foreach my $s (@seqs)
  {
    print OUT "$s\n";
  }
  close OUT;

  $fixedTxFile = $outDir . "/spp_new.tx";
  if (-e $fixedTxFile)
  {
    $TxFile = $fixedTxFile;
    if ($debug)
    {
      print "\nWorking on genus $genus w/taxonomy from $fixedTxFile\n";
    }
  }
  else
  {
    if ($debug)
    {
      print "\nWorking on genus $genus w/taxonomy from $TxFile\n";
    }
  }
  $fixedLinFile = $outDir . "/spp_new.lineage";
  if (-e $fixedLinFile)
  {
    $LineageFile = $fixedLinFile;
    if ($debug)
    {
      print "\nWorking on genus $genus w/taxonomy from $fixedLinFile\n";
    }
  }
  else
  {
    if ($debug)
    {
      print "\nWorking on genus $genus w/taxonomy from $LineageFile\n";
    }
  }
  @misclass = @{$genera{$genus}};
  $cmd = "merging.pl -t $TxFile -m $misclassFile -p $pecanFile -l $LineageFile -o $outDir $debugStr";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
}



####################################################################
##                               SUBS
####################################################################
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

exit;
