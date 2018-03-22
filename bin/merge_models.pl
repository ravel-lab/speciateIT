#!/usr/bin/env perl

=head1 NAME

  merge_models.pl

=head1 DESCRIPTION

  Using the classification scores produced by score_classification.pl, determine which PECAN 
  models can be merged by determining the reference species that are misclassified, which species the
  reference sequences are misclassified to, and merging the models together by re-labeling the
  associated taxonomy and altering the species annotation in the lineage file.

  1. For each sequence, determine if the PECAN and reference classification match. 
      If they do not match, push to hash an array with reference species as the key, and an array containing a list of species
      that the reference was misclassified to. ref -> mis1, mis2, mis3, ... , misN
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

## use Math::Permute::Array;

$OUTPUT_AUTOFLUSH = 1;
####################################################################
##                             OPTIONS
####################################################################
GetOptions(
  "output-dir|o=s"      => \my $outDir,
  "lineage-file|l=s"    => \my $LineageFile,
  "taxonomy-file|t=s"   => \my $TxFile,
  "class-scores|c=s"    => \my $clScoreFile,
  "pecan-output|p=s"    => \my $pecanFile,
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

elsif ( !$LineageFile && !$TxFile && !$clScoreFile && !$pecanFile && !$outDir )
{
  warn "\n\n\tERROR: No input, output directories or files provided.";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 0;
}

elsif ( !$outDir )
{
  warn "\n\n\tERROR: No input directory provided. Results written to merge_out";
  $outDir = "merge_out";
}

if ( !$LineageFile )
{
  warn "\n\n\tERROR: No input lineage file provided.";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 0;
}

if ( !$TxFile )
{
  warn "\n\n\tERROR: No input taxonomy file provided.";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 0;
}

if ( !$clScoreFile )
{
  $clScoreFile = "clScore.txt";
}

if ( !$pecanFile)
{
  $pecanFile = "MC_order7_results.txt";
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

####################################################################
##                               MAIN
####################################################################

print "\n---Using $LineageFile for lineage file\n";
print "---Using $TxFile for taxonomy file\n";
print "---Using $outDir for writing new files\n";

if (!-d $outDir)
{
  my $cmd = "mkdir $outDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
}

## Read in scores, taxonomy, and PECAN output file
my %lineage = read2colTbl($LineageFile);
my %clScore = readTbl($clScoreFile);
my %tx = readTbl($TxFile);
my %pecan = readTbl($pecanFile); 
my @misclass;
my $refsp;
my $pecansp;
my %newLin;
my $appended;
my %newTx;

my $misclassFile = "misclassified_taxa.txt";
open MIS, ">$misclassFile" or die "Cannot open $misclassFile for reading: $OS_ERROR\n";

foreach my $seqID (keys %clScore) ## For each of the sequences
{
  my $score = $clScore{$seqID};
  if ($score == "0")
  {
    push @misclass, $seqID;
    print MIS "$seqID\n";
  }
}
close MIS;
print "---" . scalar @misclass . " sequences written to $misclassFile\n";

merge_models(\@misclass, \%tx, \%lineage, \%pecan);


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

sub merge_models
{
  my $rseqs_to_merge = shift;
  my @seqs_to_merge = @$rseqs_to_merge;
  my $rnewtx = shift;
  my %newtx = %$rnewtx;
  my $rlineage = shift;
  my %lineage = %$rlineage;
  my $rpecan = shift;
  my %pecan = %$rpecan;
  my %doneSeqs;
  my %misTx;
  my $tx;
  my %doneTaxa;
  my $count=0;

## for all misclassified sequences, get their classified taxonomy, 
## push to a table that taxon and the sequence.
  foreach my $s (@seqs_to_merge) 
  {
    $tx = $pecan{$s};
    push @{ $misTx{$tx} }, $s;
  }
  print "---Evaluating " . scalar @seqs_to_merge ." misclassified sequences for model merging.\n\n";
  foreach my $tx ( sort { scalar(@{$misTx{$b}}) <=> scalar(@{$misTx{$a}}) } keys %misTx )
  {
    print "\nWorking on $tx (".scalar(values @{$misTx{$tx}})." sequences)\n";
    my $count = 0;
    my @clSeqs;
    my $newsp;
    my $refsp;
    my $l;
    my $clsp;
    my @actualseqID;
    my @clseqID;
    my $nrefseq;
    my $NactualSeq;
    my %cmp;
    my @refseqID;
    my @clseqs;
    my @refsps;

    $clsp = $tx;
    @clseqID = @{$misTx{$tx}};
    if ($debug)
    {
      print "There are ". scalar(@clseqID) . " sequences misclassified to $tx\n";
    }

    if ($clsp =~ /^g_/ || $clsp =~ /^f_/ || $clsp =~ /^o_/ || $clsp =~ /^c_/ || $clsp =~ /^p_/ || $clsp =~ /^d_/)
    {
      if ($debug)
      {
        print "Skipping $clsp...\n";
      }
      $count++;
      next;
    }

    #get all of the sequences that misclassified to the taxon.
    ## Occurance table
    foreach my $x (@clseqID) 
    {
      $refsp = $newtx{$x}; ## use the seqID to get its current reference species.
      chomp $refsp;
      if ($refsp =~ /$clsp/) ## If the reference species == clsp, skip this sequence
      {
        next;
      }
      else
      {
        $cmp{$refsp}++; ## If the reference species ne clsp, add refsp to cmp tbl as key, and increase value by 1.
      }
    }
    if ($debug && %cmp)
    {
      print "cmp table for $clsp: \n";
      print Dumper \%cmp;
    }
    foreach my $x (keys %cmp) ## for each ref sp that clsp belong to
    {
      $refsp = $x; 
      if ($doneTaxa{$refsp}) ## update the refsp in case it changed in a previous round.
      {
        print "refsp $refsp now $doneTaxa{$x}\n";
        $refsp = $doneTaxa{$x};
      }
      if ($doneTaxa{$clsp}) ## update the refsp in case it changed in a previous round.
       {
         print "clsp $clsp now $doneTaxa{$clsp}\n";
         $clsp = $doneTaxa{$clsp};
       }
      @clseqs = grep {$newtx{$_} =~ /$refsp/} @clseqID; ## get all of the clSeqIDs of this refsp
      
      if ($debug)
      {
        print "clseqs array (should equal n $refsp refseqs classified to $clsp: \n";
        print Dumper \@clseqs;
      }
      
      @refseqID = grep { $newtx{$_} =~ /$refsp/ } keys %newtx; ## get all of the seqIDs of the refsp
      $nrefseq = scalar(@refseqID);
      if ($debug)
      {
        print "\trefsp $refsp has $nrefseq sequence(s)\n";
      }
      if ($refsp =~ /$clsp/)
      {
        next;
      }
      elsif ($cmp{$x} <= $nrefseq*0.5) ## if the ref. sp. exists in the occurrance table
      {
        foreach my $seq (@clseqs)
        {
          $newtx{$seq} = $clsp;
          if ($debug)
          {
            print "\t$seq is now $clsp\n";
          }
        }
      }
      elsif ($cmp{$x} > $nrefseq*0.5)
      {
        my ($g1, $sp1) = split (/_/, $refsp);
        my ($g2, $sp2) = split (/_/, $clsp);
        if ($g1 =~ /$g2/)
        {
          # if ($doneTaxa{$clsp}) ## update clsp in case it has been changed to have concatenation
          # {
          #   print "\tclsp $clsp now $doneTaxa{$clsp}\n";
          #   $clsp = $doneTaxa{/$clsp/};
          # }
          $newsp = "$refsp"."_"."$clsp";
          chomp $newsp;
          $lineage{$newsp}=$lineage{$tx};
          # delete $doneTaxa{$refsp};
          # delete $doneTaxa{$clsp};
          if ($debug)
          {
            print "\t$refsp & $clsp now $newsp\n";
          }
          @actualseqID = grep { $newtx{$_} =~ /$clsp/ } keys %newtx;
          @refseqID = grep { $newtx{$_} =~ /$refsp/ } keys %newtx;
          
          ## lineages not being transferred 
          ## some things being repeated
          $doneTaxa{$clsp} = $newsp ;
          $doneTaxa{$refsp} = $newsp ;
          if ($debug)
          {
            print Dumper \%doneTaxa;
          }
          my $acount=0;
          my $rcount=0;
          my $ccount=0;
          foreach my $aseqID ( @actualseqID) ## Rename all sequences of the classified species
          {
            {
              $newtx{$aseqID} = $newsp;
              $count++;
              $acount++;
              #print "$newsp lineage becoming $l\n";
            }
          }
          foreach my $rseqID ( @refseqID) ## Rename all sequences of the reference species
          {
            {
              $newtx{$rseqID} = $newsp;
              $count++;
              $rcount++;
            }
          }
          foreach my $cseqID ( @clSeqs) ## Rename all sequences of the reference species
          {
            {
              $newtx{$cseqID} = $newsp;
              $count++;
              $ccount++;
            }
          }
        }
        else
        {
          next;
        }
      }
    }
  }
  print "Taxonomy update for $count sequences\n";

  my $fixedTxFile = $outDir . "/spp_new.tx";
  print "---Writing merged taxonomy to $fixedTxFile\n";
  open OUT, ">$fixedTxFile" or die "Cannot open $fixedTxFile for writing: $OS_ERROR";
  my %finalLin;

  my $clusterFile = $outDir . "/spp_new_clusterTbl.txt";
  print "---Writing cluster table for long taxon names\n";
  open TBL, ">$clusterFile" or die "Cannot open $clusterFile for writing: $OS_ERROR";

  my %clusters;
  my @spC;
  for my $keys (keys %newtx)
  {
    my $tx1 = $newtx{$keys};
    my $n = length($tx1);
    if ($clusters{$tx1})
    {
      my $tx2 = $clusters{$tx1};
      print OUT "$keys\t$tx2\n";
      $finalLin{$tx2} = $lineage{$tx1};
    }
    elsif ($n > 100) ## if the taxon name is > 100 characters
    { 
      my ($g, $else) = split (/_/, $tx1, 2);
      my $tx2 = $g."_Cluster"; ## shorten the name
      chomp $tx2;
      @spC = grep { $clusters{$_} =~ /$g/ } keys %clusters;
      my $l = scalar(@spC);
      $l++;
      $tx2 = $g."_Cluster_".$l;
      chomp $tx2;
      print OUT "$keys\t$tx2\n";
      $finalLin{$tx2} = $lineage{$tx1};
      $clusters{$tx1} = $tx2;
    }
    else ## if the cluster is not already taken, apply it to the seq
    {
      print OUT "$keys\t$tx1\n";
      $finalLin{$tx1} = $lineage{$tx1};
    }
  }

  foreach my $c (keys %clusters)
  {
    print TBL "$clusters{$c}\t$c\n";
  }
  close OUT;
  close TBL;

  my %finaltx = readTbl($fixedTxFile);
  my $newLineageFile = $outDir . "/spp_new.lineage";
  print "---Writing merged lineages to $newLineageFile\n";
  open OUT, ">$newLineageFile" or die "Cannot open $newLineageFile for writing: $OS_ERROR";
  my @tx;

  foreach my $x (values %finaltx)
  {
    push @tx, $x;
  }
  my @uniqTx = uniq @tx;
  
  foreach my $values (@uniqTx)
  {
    if ($finalLin{$values})
    {
      print OUT "$values\t$finalLin{$values}\n";
    }  
    elsif ($lineage{$values})
    {
      print OUT "$values\t$finalLin{$values}\n";
    }
    else
    {
      print "No lineage for $values\n";
    }
  }
  close OUT;
}


exit;