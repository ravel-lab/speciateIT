#!/usr/bin/env perl

=head1 NAME

  improve_models.pl

=head1 DESCRIPTION

  Using the classification scores produced by score_classification.pl, determine which PECAN 
  models can be improved by removing sequences with scores of 0. The script proceeds as follows: 

    1. Identify models with poorly-classified reference sequences (sub poor_classification)
        - All correct (sum of scores per taxon/nSeqs of taxon = 1)
        - Partial (0 < sum of scores per taxon/nSeqs of taxon < 1)
        - All incorrect (sum of scores per taxon/nSeqs of taxon = 0)
    2. Remove misclassified reference sequences from models if they have a score = 0.
    3. Reproduce models with above seq's removed and test their classification abilities using 
       score_classification subroutine.

  This script will produce a whole new set of models from new .tx and .lineage files.

=head1 SYNOPSIS

  improve_models.pl -i <input dir> -o <output dir>

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
  Input directory containing all.fa, spp.tx, and spp.lineage, clScore.txt, and MC_order7_results.txt files. 

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

## use Math::Permute::Array;

$OUTPUT_AUTOFLUSH = 1;
####################################################################
##                             OPTIONS
####################################################################
GetOptions(
  "input-dir|i=s"       => \my $inDir,
  "output-dir|o=s"      => \my $outDir,
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
  or pod2usage(verbose => 0,exitstatus => 1);

if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit 0;
}


if ( $inDir && !$LineageFile && !$TxFile && !$clScoreFile && !$pecanFile)
{
  $LineageFile = "$inDir/spp.lineage";
  $TxFile = "$inDir/spp.tx";
  $clScoreFile = "$inDir/clScore.txt";
  $pecanFile = "$inDir/MC_order7_results.txt";
}

elsif ( !$LineageFile && !$TxFile && !$clScoreFile && !$pecanFile && !$inDir && !$outDir )
{
  warn "\n\n\tERROR: No input, output directories or files provided.";
  exit 0;
}
elsif ( !$inDir )
{
  warn "\n\n\tERROR: No input directory provided.";
  $outDir = $inDir;
}
elsif ( !$outDir )
{
  warn "\n\n\tERROR: No input directory provided.";
  $outDir = "merge_model_out";

}

if ( $inDir && !$LineageFile )
{
  $LineageFile = "$inDir/spp.lineage";
}

if ( $inDir && !$TxFile )
{
  $TxFile = "$inDir/spp.tx";
}

if ( $inDir && !$clScoreFile )
{
  $clScoreFile = "$inDir/clScore.txt";
}

if ( $inDir && !$pecanFile)
{
  $pecanFile = "$inDir/MC_order7_results.txt";
}

my $quietStr = "";
if ( $quiet  )
{
  $quietStr = "--quiet";
}

my $debugStr = "";
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

if (!$offsetCoef)
{
  $offsetCoef = 0.9;
}

if (!$txSizeThld)
{
  $txSizeThld = 10;
}
####################################################################
##                               MAIN
####################################################################

print "\n---Using $LineageFile for lineage file\n";
print "---Using $TxFile for taxonomy file\n";
print "---Using $clScoreFile for classification score file\n";
print "---Using $outDir for writing new files\n";

if ($outDir ne $inDir)
{
  my $cmd = "rm -rf $outDir; mkdir $outDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
}

## Read in scores, taxonomy, and PECAN output file
my %clscore = readTbl($clScoreFile);
my %tx = readTbl($TxFile);
my %pecan = readTbl($pecanFile); 

## Make a unique array of taxa
my @uniqtx = uniq (values %tx);

print "---Determining models with poorly-classified sequences\n";
my ($rmisclassified, $rpoorclass) = poor_classification(\@uniqtx);

## Push to an single array all of the seq's that were misclassified (score = 0)
## For removal from dataset.

push (@$rmisclassified, @$rpoorclass);

my $rmFile = "$outDir/seqs_to_remove.seqID";
print "---Writing IDs of misclassified sequences contributing < 50% to their species to $rmFile\n";
open OUT, ">$rmFile" or die "Cannot open $rmFile for writing: $OS_ERROR";
print OUT @$rmisclassified;
close OUT;

# my $grFile = "$outDir/seqs_greater_than_50.seqID";
# print "---Writing IDs of misclassified sequences contributing > 50% to their species to $grFile\n";
# open OUT, ">$grFile" or die "Cannot open $grFile for appending: $OS_ERROR";
# print OUT @$seqs_greater_than_50;
# close OUT;

# my $eqFile = "$outDir/seqs_equal_to_50.seqID";
# print "---Writing IDs of misclassified sequences contributing = 50% to their species to $eqFile\n";
# open OUT, ">$eqFile" or die "Cannot open $eqFile for appending: $OS_ERROR";
# print OUT @$seqs_equal_to_50;
# close OUT;

print "---Removing misclassified sequences (classification score = 0) from taxonomy file\n";
my $newtxFile = "$outDir/spp_new.tx";
my $cmd = "select_tx.pl -i $TxFile -e $rmFile -o $newtxFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

print "---Removing misclassified sequences (classification score = 0) from lineage file\n";
my $newLinFile = "$outDir/spp_new.lineage";
$cmd = "select_fullTx.pl -t $newtxFile -f $LineageFile -o $newLinFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

print "--- Building new fasta file with misclassified sequences removed\n";
my $faFile = $inDir . "/all.fa";
my $reducefaFile = $outDir . "/select.fa";
$cmd = "select_seqs.pl -i $faFile -s $newtxFile -o $reducefaFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

print "--- Building model tree and creating taxon's reference fasta files\n";
$cmd = "buildModelTree $quietStr -l $newLinFile -i $reducefaFile -t $newtxFile -o $outDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

print "--- Building MC models\n";
$cmd = "buildMC -t $outDir/spp_paths.txt -k 8 -d $outDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

print "--- Estimating error thresholds\n";
$cmd = "est_error_thlds -v --offset-coef $offsetCoef --tx-size-thld $txSizeThld -d $outDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

print "--- Running classify on the ref fasta file\n";
$cmd = "classify $skipErrThldStr -d $outDir -i $reducefaFile -o $outDir"; # $ppEmbeddingStr
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

print "--- Comparing ref seq's taxonomy with the classification results\n";
my $sppSeqID = "$outDir/sppSeqID.lineage";
$cmd = "get_lineage.pl -a $newtxFile -l $newLinFile -o $sppSeqID";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;


print "---Using $sppSeqID for lineage file\n";
my %newlineage = read2colTbl($sppSeqID);

my $newpecanFile = "$outDir/MC_order7_results.txt";
my %newpecan = readTbl($newpecanFile); 
print "---Using $newpecanFile for lineage file\n";

my @clIDs = keys %newpecan;

print "--- Scoring ref seq's taxonomy with the new classification results\n";
my $rclScore = score_classification( \@clIDs, \%newlineage, \%newpecan);
my %clScore = %$rclScore;
my $clScore_sum = sum values %clScore;
my $clScore_perc = 100.0 * $clScore_sum /  scalar @clIDs; 

print "\n\nUsing error thresholds, the new classification score is: " . $clScore_perc . "\n";
my $newclScoreFile = "$outDir/clScore_new.txt";
open OUT, ">$newclScoreFile" or die "Cannot open $newclScoreFile for appending: $OS_ERROR";
print OUT "seqID\tclassification_score\n";
foreach my $key (keys %clScore)
{
  print OUT $key . "\t". $clScore{$key} . "\n";
}
close OUT;

exit;

# print "---Using $newtxFile for taxonomy file\n";
# my %newtx = readTbl($newtxFile);

# print "---Using $newLinFile for lineage file\n";
# my %lineage = read2colTbl($newLinFile);
# my @poorclass = @$rpoorclass;

# print "---Merging species of misclassified sequences which contributed to >= 50% of their species\n";
# my ($rfixedtx, $rnewlineage, $rdeletedseqs) = merge_models(\%newtx, \%lineage, \%pecan, \@poorclass);

# # my %fixedtx = %$rfixedtx;
# # my %newlineage = %$rnewlineage;
# # my @deletedSeqs = @$rdeletedseqs;

# my $fixedTxFile = "$outDir/spp_merged.tx";
# print "---Writing merged taxonomy to $fixedTxFile\n";
# open OUT, ">$fixedTxFile" or die "Cannot open $fixedTxFile for appending: $OS_ERROR";
# for my $keys (keys %fixedtx)
# {
#   print OUT "$keys\t$fixedtx{$keys}\n";
# }
# close OUT;

# print "--- Using $fixedTxFile to obtain lineage file\n";
# my $newLineageFile = "$outDir/spp_merged.lineage";
# $cmd = "select_fullTx.pl -t $fixedTxFile -f $LineageFile -o $newLineageFile";
# print "\tcmd=$cmd\n" if $dryRun || $debug;
# system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
# print "---Writing merged lineages to $newLineageFile\n";
# open OUT, ">$newLineageFile" or die "Cannot open $newLineageFile for appending: $OS_ERROR";
# for my $keys (keys %newlineage)
# {
#   print OUT "$keys\t$newlineage{$keys}\n";
#   if (!$debug)
#   {
#     print "Printed to $newLineageFile: \t $keys\t$newlineage{$keys}\n";
#   }
# }
# close OUT;

print "---Removed ". scalar @$rmisclassified . " misclassified sequences from taxonomy & lineage files because they were < 50% of their ref spp.\n";
#print "---Proceed to rebuilding models using $fixedTxFile and $newLineageFile as input\n";
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

sub poor_classification{
  my $runiqtx = shift;
  my @uniqtx = @$runiqtx;
  my $spp=0;
  my @seqIDs;
  my %sppSeqID;
  my %spsums;
  my $score;
  my $sum;
  my $avg;
  my $id;
  my @seqs_to_remove;
  my @seqs_greater_than_50;
  my @seqs_equal_to_50;
  my @misclassified;
  my @poorclass;


  ## For each unique taxon of the reference data 
  ## Obtain all seqIDs of that taxon...
  ## and build array of corresponding classification scores
  foreach my $key ( @uniqtx) 
  {
    @seqIDs = grep { $tx{$_} eq $key } keys %tx; 
    my @clscores;
    foreach $id ( @seqIDs ) 
    {
      if ( exists $clscore{$id} ) 
      {
        $score = $clscore{$id}; 
        push @clscores, $score; 
        if ($debug)
        {
          print "The score for $id was $clscore{$id}\n";
        }
      }
    }

    ##Place the clScore sum into a hash for that species.
    $sum = sum @clscores;
    $avg = $sum / scalar @clscores;
    $spsums{$key} = $avg;
    $spp++;

    ## If the mean classification score for a taxon is eq 1, 
    ## Go to the next taxon
    if ( $avg == 1 )
    {
      if ($debug)
      {
        print "All sequences of $key are correctly classified\n";
      }
    }
    ## If the average classification score of a taxon is less than 1
    ## loop through the sequences of that taxon and if their score
    ## is 0, push the seqID to the misclassified array
    elsif ( $avg < 1 )
    {
      foreach $id ( @seqIDs ) 
      {
        $score = $clscore{$id};
        if ($score == 0)
          {
          push @misclassified, $id."\n";
          }
        else
        {
          push @poorclass, $id."\n";
        }
      }
      #print "Only some of the sequence(s) (" . scalar @clscores . " total) of $key are misclassified\n";
      #print scalar @misclassified . " sequences have a score of 0\n"; 

    ## Determine the number of misclassified sequences in a reference taxon
    ## If the the difference is < 50%, push the misclassified seq IDs to 
    ## an array for removal
    # my $diff = scalar @misclassified / scalar @seqIDs;
    # if ($diff < 0.5)
    # {
    #   foreach my $x (@misclassified)
    #   {
    #     push @seqs_to_remove, $x."\n";
    #   }
    # }
    # ## If the the difference is > 50%, push the misclassified seq IDs to 
    # ## an array
    # elsif ($diff > 0.5)
    # {
    #   foreach my $x (@misclassified)
    #   {
    #     push @seqs_greater_than_50, $x."\n";
    #   }
    # }
    # elsif ($diff = 0.5)
    # {
    #   foreach my $x (@misclassified)
    #   {
    #     push @seqs_equal_to_50, $x."\n";
    #   }
    # }
  }
}
return (\@misclassified, \@poorclass)
}

sub merge_models{
 my ($rnewtx, $rlineage, $rpecan, $rseqs_to_merge) = @_;
  
  my %newtx = %{$rnewtx};
  my %lineage = %{$rlineage};
  my %pecan = %{$rpecan};
  my %doneSeqs;
  my @deletedSeqs;

  
  #my $merged_spp = "$outDir/merged_spp.log";
  #open OUT, ">$merged_spp" or die "Cannot open $merged_spp for appending: $OS_ERROR";
  #print OUT "newSpp\trefSpp\tpecanspp\n";
  foreach my $seq (@$rseqs_to_merge)
  {
    chomp $seq;
    my @clSeqs;
    my $winner;
    my $newsp;
    my $refsp;
    my $pecansp;
    my @refspp;
    my @pecansp;
    my @refseqID;
    my @misclseq;
    my $nrefseq;
    my $nclseq;
    my $r;
    my $c;
    my @urefsp;
    my $splin;
    my $i;
    my $a;
 
    if (!exists $doneSeqs{$seq})
    {
      $pecansp = $pecan{$seq};
      {    
        print "\n\nThe species of $seq is $newtx{$seq} but was classified to $pecansp.\n";
      }
   
      ## For each sequence, identify the species it was classified as and the species it came from
      ## compare the number of sequences for each of these species
      ## and which ever one has more sequences, change the lineage (ph - sg) to the one with more seq's
      ## and change the species name for all to contenate as biggersp_smallersp
      ## Therefore, we need the species it was classified as (pecanFile), the reference (%tx), the number
      ## of seq's in each of those species (from the ref %tx)

      
      if (!$pecansp)
      {
        print "$seq not found in PECAN output.\n";
        exit;
      }

      ## Get all seqIDs from classified taxa that equal the classified species of interest in this loop
      @clSeqs = grep { $pecan{$_} eq $pecansp } keys %pecan;
      if ($debug)
      {
        print "There are " . scalar @clSeqs . " sequences in MC_output for $pecansp\n";
      }
      
      ## Now use those seqIDs to build an array of all of their reference species
      foreach my $clseq ( @clSeqs)
      {
        chomp $clseq;
        my $a = $newtx{$clseq};
        push @refspp, $a;
      }

      if (@refspp)
      {
        @urefsp = uniq @refspp;
      }
      else 
      {
        print "array refspp for $pecansp is empty\n";
      }

        ## If the reference sp has more seq's, then we need to 
        ## 1. update the taxonomy (seqid -> tx) to be refsp_pecansp
        ## 2. delete the pecansp from the lineage hash
        ## 3. add a new key/value with the refsp_pecansp to the ref_lineage

      ## If there is only one reference species for the seq's that misclassified, merge the two species (all seq's of both species become same taxon)
      if (@urefsp && scalar @urefsp == 2)  
      {
        $refsp = $newtx{$seq};
        @refseqID = grep { $newtx{$_} eq $refsp } keys %newtx; ## get all seqIDs of this reference species.
        if ($debug)
          {
          print "There are " . scalar @refseqID . " sequences in $refsp\n";
          }
        $nrefseq = scalar @refseqID;

        ## Determine how many of the seq's from the reference species are misclassified. If > 2 classifications, remove reference species (seqs and taxon)
        my @rseqs_miscl;
        foreach my $r (@refseqID)
        {
          $c = $pecan{$r};
          push @rseqs_miscl, $c; ## push to the array the classified species of that sequence
        }
        my @urseqs_miscl = uniq @rseqs_miscl; ## produce a unique version of this array
        if (scalar @urseqs_miscl > 2) ## if the array has more than 2 species (1 wrong and 1 right), then remove the reference species entirely.
        {
          foreach my $r (@refseqID)
          {
           delete $newtx{$r};
          }
          push @deletedSeqs, @refseqID;
          delete $lineage{$refsp};
          print "Sequences of $refsp were misclassified to multiple species. Removed $refsp seq's from taxonomy, and $refsp from lineage.\n";
        }
        if (scalar @refseqID < 10)
        {
          foreach my $r (@refseqID)
          {
           delete $newtx{$r};
          }
          push @deletedSeqs, @refseqID;
          delete $lineage{$refsp};
          print "$refsp contained only ". scalar @refseqID . " sequences. Removed $refsp seq's from taxonomy, and $refsp from lineage.\n";
        }

        %doneSeqs = map { $_ => 1 } @refseqID;  

        # @misclseq = grep { $newtx{$_} eq $pecansp } keys %newtx; ## get all seqIDs of the classified species
        # if ($debug)
        #   {
        #   print "There are " . scalar @misclseq . " sequences in $pecansp\n";
        #   }
        # $nclseq = scalar @misclseq;

        # $r = scalar @refseqID;
        # $c = scalar @misclseq;

        # if ( $r > $c)
        # {
        #   $winner = 0;
        # }
        # elsif ( $r < $c)
        # {
        #   $winner = 1;
        # }
        # elsif ( $r eq $c)
        # {
        #   $winner = 2;
        #   print "***WARNING: There are the same number of seq's in $refsp ($r) and $pecansp ($c)\n";
        # }
        # ## If there are more seq's in the reference species than the classified species, the newSp starts with the ref name
        # if ($winner == 0)
        # {
        #   ## define the new spp merged name
          
        #   foreach my $clseq ( @misclseq)
        #   {
        #     chomp $clseq;
        #     delete $newtx{$clseq};
        #     delete $newtx{$seq};
        #   }          

        #   delete $lineage{$pecansp};
        #   #delete $lineage{$refsp};
        #   #$lineage{$newsp} .= $splin;

        #   if ($debug)
        #   {
        #     print "$pecansp deleted from the lineage file, and contained $c sequences which were deleted from the taxonomy file.\n";
        #   }

        # }
        # elsif ($winner >= 1)
        # {
        #   $newsp = "$pecansp"."_"."$refsp";
        #   chomp $newsp;
        #   my @t = split "_", $newsp;

        #   if ( @t > 3 && $t[$#t] ne "etal" )
        #   {
        #     $newsp = join "_", @t[0..2];
        #     $newsp .= "_etal";
        #   }

        #   $splin = $lineage{$pecansp};
        #   delete $lineage{$pecansp};
        #   #delete $lineage{$refsp};
        #   $lineage{$newsp} .= $splin;
        #   delete $newtx{$seq};
        #   foreach my $clseq ( @misclseq)
        #   {
        #     chomp $clseq;
        #     $newtx{$clseq} = $newsp;
        #   }

        #   print "$pecansp is now $newsp. Name changed for " . scalar @misclseq . " sequences.\n";
        # }
       }
       elsif (@urefsp && scalar @urefsp < 2)
       {
        print "As per $seq, sequences classified to $pecansp come from fewer than 2 reference species (" . scalar @urefsp . "). Investigate this.\n";
       }
       else ## If there are multiple reference species all classifying to the same species (incorrectly), remove those sequences only.
        {
          print "Sequences classified to $pecansp come from " . scalar @urefsp . " reference species. Deleting $seq. \n";
          delete $newtx{$seq};
          push @deletedSeqs, $seq;
          $doneSeqs{$seq} = 1 ;   
        }
      #print OUT "$newsp\t$refsp\t$pecansp\n";
      push (@misclseq, @refseqID);
      push (@pecansp, @refspp);
        ## For the sequences classified to only this species coming from only one species, we need to change the taxonomy for
        ## both the clSeqs and the refSeqs, produce a lineage line for the newspecies, and remove from lineage the pecansp and the refSp
      
      #%doneSeqs = map { $_ => 1 } @clSeqs;
      #print "Sequences misclassified to $pecansp completed\n\n";
    }
  }
  close OUT;
  #print "Species changes captured in $merged_spp\n";
  return (\%newtx, \%lineage, \@deletedSeqs);
}

sub score_classification{ 
  my $rclIDs = shift;
  my $rtestLineage = shift;
  my $rclTx = shift;

  my %testLineage = %$rtestLineage;
  my %clTx = %$rclTx;

  my %clScore;

  my @clIDs = @$rclIDs;
  foreach my $s (@clIDs) ## For each classified sequence
  {
    chomp $s;
    ##print "\nHere is the lineage for $s: ". $testLineage{$s} . "\n";

    ## t is equal to the classified taxa
    my $t = $clTx{$s};
    if (!$t)
    {
      next;
    } 
    else
    {
      ##print "The value of $s is $clTx{$s}.\n";
    }
    ## a is equal to the reference lineage
    my $a = $testLineage{$s};
    if (!$a)
    {
      if ($debug)
      {
        print "\nWARNING: $s not found in testLineage hash table\n;"
      }
      next;
    }
    #print "\nThe value of $s is $testLineage{$s}\n";
    ## if the classified taxa is found in the reference lineage
    my $value = 0;
    if ( grep /$t/, $a) 
    {   
      #print "\n $t is part of $a.\n";
      $value++;
      if ($debug)
      {
        print "\n $t is part of $a.\n";
      }
      if ( $value > 0) ##If the classifed taxonomy is found in the array
      {
        if ($clTx{$s} =~ /^d_/) 
        {
          $clScore{$s} = 0;
        }
        elsif ($clTx{$s} =~ /^p_/) 
        {
          $clScore{$s} = 0.25;
        }
        elsif ($clTx{$s} =~ /^c_/) 
        {
          $clScore{$s} = 0.325;
        }
        elsif ($clTx{$s} =~ /^o_/) 
        {
          $clScore{$s} = 0.5;
        }
        elsif ($clTx{$s} =~ /^f_/) 
        {
          $clScore{$s} = 0.625;
        }
        elsif ($clTx{$s} =~ /^g_/) 
        {
          $clScore{$s} = 0.75;
        }          
        elsif ($clTx{$s} =~ /^sg_/) 
        {
          $clScore{$s} = 0.825;
        }
        else 
        {
          $clScore{$s} = 1;
        }
      }
      elsif ($value = 0) ## And of course, these have to be correct with the lineage. If they are not correct, they get a zero.
      {
        $clScore{$s} = 0;
      }   
      #print "\nThe score for $s is $clScore{$s}\n";
    }
    else ## And of course, these have to be correct with the lineage. If they are not correct, they get a zero.
    {
      $clScore{$s} = 0;
    } 
  }
return \%clScore;
}

exit;