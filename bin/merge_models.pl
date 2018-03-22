#!/usr/bin/env perl

=head1 NAME

  merge_models.pl

=head1 DESCRIPTION

  Using the classification scores produced by score_classification.pl, determine which PECAN 
  models can be merged by determining the reference species that are misclassified, which species the
  reference sequences are misclassified to, and merging the models together by re-labeling the
  associated taxonomy and altering the species annotation in the lineage file.

  1. For each sequence, determine if the PECAN and reference classification match. 
      -Use clScore.txt to get seqIDs w/score of 0
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

####################################################################
##                               MAIN
####################################################################
my $outDir = "merge_out";
print "---Using $LineageFile for lineage file\n";
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

  my $fixedTxFile = $outDir . "/spp_new.tx";
  if (-e $fixedTxFile)
  {
    unlink $fixedTxFile;
  }

## for all misclassified sequences, get their classified taxonomy, 
## push to a table that taxon and the sequence.
  foreach my $s (@seqs_to_merge) 
  {
    $tx = $pecan{$s};
    push @{ $misTx{$tx} }, $s;  ## %misTx: classified taxon => array of seqIDs
  }
  print "---Evaluating " . scalar @seqs_to_merge ." misclassified sequences for model merging.\n\n";
  foreach my $tx ( sort { scalar(@{$misTx{$b}}) <=> scalar(@{$misTx{$a}}) } keys %misTx ) ## sort %misTx by the number of misclassified seq's 
  {
    print "\nWorking on $tx\n"; ## (".scalar(values @{$misTx{$tx}})." sequences misclassified to $tx)\n";
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
   
    print "\tThere were ". scalar(@clseqID) . " sequences misclassified to $tx\n";
    
    ## Evaluate if clsp is anything but a species. 
    if ($clsp =~ /^g_/ || $clsp =~ /^f_/ || $clsp =~ /^o_/ || $clsp =~ /^c_/ || $clsp =~ /^p_/ || $clsp =~ /^d_/)
    {
      if ($debug)
      {
        print "\tSkipping $clsp...\n";
      }
      $count++;
      next;
    }

    ## Occurance table %cmp: $refsp => nSeqs classified to clSp that are refsp
    foreach my $x (@clseqID) 
    {
      $refsp = $newtx{$x}; ## use the seqID to get its current reference species.
      chomp $refsp;
      if ($refsp =~ /$clsp/) ## If the reference species == clsp, skip this sequence SHOULD NEVER HAPPEN
      {
        print "encountered refsp == clsp -- CHECK\n";
        next;
      }
      else
      {
        $cmp{$refsp}++; ##  %cmp: $refsp => nSeqs classified to clSp that are refsp
      }
    }
    if ($debug && %cmp)
    {
      print "Hash $clsp cmp:\n";
      print Dumper \%cmp;
    }
    foreach my $x (keys %cmp) ## for each ref sp that clsp belong to
    {
      ### UPDATE TAXONOMY IF CHANGED DUE TO MODEL MERGING ##### 
      if (-e $fixedTxFile)
      {
        print "---Writing merged taxonomy to $fixedTxFile\n";
        open OUT, "<$fixedTxFile" or die "Cannot open $fixedTxFile for reading: $OS_ERROR";
        my %newtx = readTbl($fixedTxFile);
        close OUT;
      }
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
      if ($refsp =~ /$clsp/ || $clsp =~ /$refsp/)
      {
        next;
      }

      ### GET ARRAYS OF SEQIDS FOR REFSP AND ALL REFSP SEQ'S MISCLASSIFIED TO CLSP ##### 
      @clseqs = grep {$newtx{$_} =~ /$refsp/} @clseqID; ## for this refsp, get all seqIDs that misclassified to clSp
      if ($debug)
      {
        print "Array $refsp clSeqs:\n";
        print Dumper \@clseqs;
      }
      
      @refseqID = grep { $newtx{$_} =~ /$refsp/ } keys %newtx; ## get all of the seqIDs of this refsp from the reference taxonomy.
      $nrefseq = scalar(@refseqID);
      if ($debug)
      {
        print "\trefsp $refsp has $nrefseq sequence(s)\n";
      }

      ### COMPARE ARRAYS OF SEQIDS FOR REFSP AND ALL REFSP SEQ'S MISCLASSIFIED TO CLSP ##### 
      if ($cmp{$x} <= $nrefseq*0.5) ## If n RefSeqs classified to clSp that are <= 50% of nRefSeqs 
      {
        foreach my $seq (@clseqs)
        {
          $newtx{$seq} = $clsp; ## Just rename those misclassified seq's to the clSP
          print "\t$seq is now $clsp\n";
        }
      }
      elsif ($cmp{$x} > $nrefseq*0.5) ## If n RefSeqs classified to clSp that are > 50% of nRefSeqs ... merge if genera are the same
      {
        my ($g1, $sp1) = split (/_/, $refsp);
        my ($g2, $sp2) = split (/_/, $clsp);
        if ($g1 =~ /$g2/)
        {
          $newsp = "$refsp"."_"."$clsp"; ## Concatenate names
          chomp $newsp;
          $lineage{$newsp}=$lineage{$tx}; ## Produce lineage for concatenation
          print "\t$refsp & $clsp => $newsp\n";

          ### KEEP RECORD OF OLD NAMES AND NEW NAME
          $doneTaxa{$clsp} = $newsp ;
          $doneTaxa{$refsp} = $newsp ;
          if ($debug)
          {
            print "Hash doneTaxa:\n";
            print Dumper \%doneTaxa;
          }

          ### RENAME ALL SEQS OF EACH SPECIES TO THE NEWLY CONCATENATED NAME
          @actualseqID = grep { $newtx{$_} =~ /$clsp/ } keys %newtx;
          @refseqID = grep { $newtx{$_} =~ /$refsp/ } keys %newtx;
          
          my $count=0;
          foreach my $aseqID ( @actualseqID) ## Rename all sequences of the classified species
          {
            {
              $newtx{$aseqID} = $newsp;
              $count++;
            }
          }
          foreach my $rseqID ( @refseqID) ## Rename all sequences of the reference species
          {
            {
              $newtx{$rseqID} = $newsp;
              $count++;
            }
          }
        }
        else
        {
          print "Genera do not match... skipping\n";
          next;
        }
      }
    }
    print "---Writing merged taxonomy to $fixedTxFile\n";
    if (-e $fixedTxFile)
    {
      unlink $fixedTxFile;
    }
    open OUT, ">$fixedTxFile" or die "Cannot open $fixedTxFile for writing: $OS_ERROR";
    for my $key (keys %newtx)
    {
      print OUT "$key\t$newtx{$key}\n";
    }
    close OUT;
  }
  print "Taxonomy updated for $count sequences\n";

  print "---Writing final merged taxonomy to $fixedTxFile\n";
  if (-e $fixedTxFile)
  {
    unlink $fixedTxFile;
  }
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