#!/usr/bin/env perl

=head1 NAME

  merging.pl

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

    merging.pl -t $TxFile -m $misclassFile -p $pecanFile -l $LineageFile

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
  "misclass-file|m=s"   => \my $misclassFile,
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

if ( !$LineageFile && !$TxFile && !$misclassFile && !$pecanFile)
{
  print "Please provide lineage (-l), taxonomy (-t), misclass (-m), and PECAN (-p) files\n";
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
elsif ( !$misclassFile )
{
  print "Please provide genus misclassified (-t) files\n";
  exit 0;
}

elsif ( !$pecanFile)
{
  $pecanFile = "MC_order7_results.txt";
}
if (!$outDir)
{
  my $outDir = "merge_out";
}
#print "Perl $^V\n";

####################################################################
##                               MAIN
####################################################################
# print "\tUsing $LineageFile for lineage file\n";
# print "\tUsing $TxFile for taxonomy file\n";
# print "\tUsing $outDir for writing new files\n";

## Read in scores, taxonomy, and PECAN output file
my %lineage = read2colTbl($LineageFile);
my %pecan = readTbl($pecanFile);
my @misclass = readArray($misclassFile);
my %newtx = readTbl($TxFile);
my $fixedTxFile = $outDir . "/spp_new.tx";

if (-e $fixedTxFile)
{
  %newtx = readTbl($fixedTxFile);
}
else 
{
  %newtx = readTbl($TxFile);
}

my $clusterFile = $outDir . "/spp_new_clusterTbl.txt";
my %clusters;
if (-e $clusterFile)
{
  %clusters = readTbl($clusterFile);
}

## for all misclassified sequences, get their classified taxonomy, 
my %misTx;
my $tx;
foreach my $s (@misclass) 
{
  $tx = $pecan{$s};
  push @{ $misTx{$tx} }, $s;  ## %misTx: classified taxon => array of seqIDs
}

my %donetaxa;

### PERFORM COMPARISONS AND ANNOTATION UPDATES
foreach my $x (sort keys %misTx)
{
  if ($debug)
  {
    print "Working on $x\n"; ## (".scalar(values @{$misTx{$tx}})." sequences misclassified to $tx)\n";
  }
  my $newsp;
  my $l;
  my @actualseqID;
  my @clseqID;
  my $nrefseq;
  my $NactualSeq;
  my %cmp;
  my @refseqID;
  my @clseqs;
  my $refsp;
  my $sp;

  @clseqID = @{$misTx{$x}};
  if ($debug)
  {
    print "\tThere were ". scalar(@clseqID) . " sequences misclassified to $x\n";
  }
  ## Evaluate if clsp is anything but a species. 
  
  ## Occurance table %cmp: $refsp => nSeqs classified to clSp that are refsp
  foreach my $y (@clseqID) 
  {
    $refsp = $newtx{$y}; ## use the seqID to get its current reference species.
    chomp $refsp;
    if ($refsp !~ /$x/ && $x !~ /^g_/ && $x !~ /^f_/ && $x !~ /^o_/ && $x !~  /^c_/ && $x !~  /^p_/ && $x !~ /d_Bacteria/)
    {
      $cmp{$refsp}++; ##  %cmp: $refsp => nSeqs classified to clSp that are refsp
    }
  }

  ## for each ref sp that clsp belong to determine how to change annotation, if needed
  foreach my $y (sort keys %cmp) 
  {
    my $clsp = $x;
    my $refsp = $y;
    my $rcluster;
    my $cluster;
    if ($debug)
    {
      print "\n\tComparing RefSp: $refsp & ClSp: $clsp\n";
    }
    ### UPDATE ANNOTATIONS IF THEY HAVE BEEN ALTERED (SHOULD BE PRESENT IN DONETAXA) ##### 
    if ($donetaxa{$refsp}) ## update the refsp in case it changed in a previous round. %donetaxa: old name => new name 
    {
      my $n = $donetaxa{$refsp};
      $refsp = $newtx{$n};
      #delete donetaxa{$sp};
      if ($debug)
      { 
        print "\t\tRefSp now $refsp\n";
      }
    }
    if ($refsp =~ /Cluster/)
    {
      $rcluster = $refsp;
      if (%clusters)
      {
        $refsp = $clusters{$rcluster};
        if ($debug)
        {
          print "\t\tRefSp now $refsp\n";
        }
      }
      else
      {
        print "\n\tERROR: Taxon is cluster, but no cluster file provided\n";
        exit;
      }
      
    }
    if ($donetaxa{$clsp}) ## update the refsp in case it changed in a previous round.
    {
       my $n = $donetaxa{$clsp};
       $clsp = $newtx{$n};
       #delete donetaxa{$y};
       if ($debug)
        {
          print "\t\tClSp now $clsp\n";
        }
    }
    if ($clsp =~ /Cluster/)
    {
      $cluster = $clsp;
      $clsp = $clusters{$cluster};
      if ($debug)
        {
          print "\t\tClSp now $clsp\n";
        }
    }

    ### GET ARRAYS OF SEQIDS FOR REFSP AND ALL REFSP SEQ'S MISCLASSIFIED TO CLSP ##### 
    @clseqs = grep {$newtx{$_} =~ /$refsp/} @clseqID;
    @refseqID = grep { $newtx{$_} =~ /$refsp/ } keys %newtx; ## get all of the seqIDs of this refsp from the reference taxonomy.
    $nrefseq = scalar(@refseqID);
    if ($debug)
    {
      print "\t\tRefSp $refsp has $nrefseq sequence(s) and " . scalar(@clseqs) . " were misclassified to $clsp\n";
    }

    ### COMPARE ARRAYS OF SEQIDS FOR REFSP AND ALL REFSP SEQ'S MISCLASSIFIED TO CLSP ##### 
    if (scalar(@clseqs) == 0)
    {
      next;
    }
    elsif ($cmp{$y} <= $nrefseq*0.05) ## If n RefSeqs classified to clSp that are <= 50% of nRefSeqs 
    {
      foreach my $seq (@clseqs)
      {
        if (defined $cluster)
        {
          $clsp = $cluster;
        }
        $newtx{$seq} = $clsp;
      }
      if ($debug)
        {
          print "\t\t" . scalar @clseqs . " sequence(s) is/are now $clsp\n";
        }
    }
    elsif ($cmp{$y} > $nrefseq*0.05) ## If n RefSeqs classified to clSp that are > 50% of nRefSeqs ... merge if genera are the same
    {
      my ($g1, $sp1) = split (/_/, $refsp);
      my ($g2, $sp2) = split (/_/, $clsp);
      if ($g1 =~ /$g2/)
      {
        my $count = 0;
        $newsp = "$refsp"."_"."$clsp"; ## Concatenate names
        chomp $newsp;
        if ($debug)
        {
          print "\t\t$refsp & $clsp => $newsp\n";
        }

        ### BUILD CLUSTERS IF NECESSARY ### 
        my $n = length($newsp);
        my @spC;
        if (defined $cluster)
        {
          $clusters{$cluster} = $newsp;
          if ($debug)
          {
            print "\t\tUpdated $cluster as $newsp\n";
          }
          $newsp = $cluster;
        }
        elsif (!defined $cluster && $n > 100)
        {
          my ($g, $else) = split (/_/, $clsp, 2);
          my $tx2 = $g."_Cluster"; ## shorten the name
          chomp $tx2;
          @spC = grep { $clusters{$_} =~ /$g/ } keys %clusters;
          my $l = scalar(@spC);
          $l++;
          $tx2 = $g."_Cluster_".$l;
          chomp $tx2;
          $clusters{$tx2} = $newsp;
          $newsp = $tx2;
        }

        ## ADD TO LINEAGE THE NEWSP AND LINEAGE FROM CLSP 
        ## ALWAYS ADD TO LINEAGE, DON'T DELETE, BECAUSE A FINAL LINEAGE WILL BE BUILT AT THE END 
        ## FROM THE COMPLETED NEWTX.
        $lineage{$newsp}=$lineage{$x};

        ### RENAME ALL SEQS OF EACH SPECIES TO THE NEWLY CONCATENATED NAME
        ### ADD TO HASH DONE TAXA ONLY WHEN ANNOTATIONS ARE MANIPULATED
        @refseqID = grep { $newtx{$_} =~ /$refsp/ } keys %newtx;
        if ($cluster)
        {
          @actualseqID = grep { $newtx{$_} =~ /$cluster/ } keys %newtx;
        }
        else
        {
          @actualseqID = grep { $newtx{$_} =~ /$clsp/ } keys %newtx;
        }
          
        foreach my $aseqID ( @actualseqID) ## Rename all sequences of the classified species
        {
          $newtx{$aseqID} = $newsp;
          $donetaxa{$clsp} = $aseqID;
          $count++;
        }
        foreach my $rseqID ( @refseqID) ## Rename all sequences of the reference species
        {
          $newtx{$rseqID} = $newsp;
          $donetaxa{$refsp} = $rseqID;
          $count++;
        }
        if ($debug)
        {
          print "\t\tHash donetaxa:\n\t\t";
          print Dumper \%donetaxa;
          print "\t\t$count sequences altered to $newsp\n";
        }
      }
      else
      {
        if ($debug)
        {
          print "\t\tGenera do not match. Skipping...\n";
        }
        next;
      } 
     if ($debug)
      {
        print "\tLength of donetaxa hash: ". scalar(keys %donetaxa) . "\n";
        print "\tSize of donetaxa hash: ". total_size(\%donetaxa) . "\n";
        print "\tLength of newtx hash: ". scalar(keys %newtx) . "\n";
        print "\tSize of newtx hash: ". total_size(\%newtx) . "\n";
        print "\tLength of clusters hash: ". scalar(keys %clusters) . "\n";
        print "\tSize of clusters hash: ". total_size(\%clusters) . "\n\t";
        print Dumper \%clusters;
      }
    }
  }
}

### Finalize files with all changes ####
#print "\tWriting final merged taxonomy to $fixedTxFile\n";
open OUT, ">$fixedTxFile" or die "Cannot open $fixedTxFile for writing: $OS_ERROR";
foreach my $s (keys %newtx)
{
  my $tx = $newtx{$s};
  print OUT "$s\t$newtx{$s}\n";
}
close OUT;

### Write clusters file for annotations > 100 characters long ####
if (%clusters)
{
  my $clusterFile = $outDir . "/spp_new_clusterTbl.txt";
  print "\tWriting cluster table for long taxon names\n";
  open TBL, ">$clusterFile" or die "Cannot open $clusterFile for writing: $OS_ERROR";
  foreach my $c (keys %clusters)
  {
    print TBL "$c\t$clusters{$c}\n";
  }
  close TBL;
}

### Produce lineage file with all updated taxonomy (including clusters) ####
my %finaltx = readTbl($fixedTxFile);
my $newLineageFile = $outDir . "/spp_new.lineage";
#print "\tWriting merged lineages to $newLineageFile\n";
open OUT, ">$newLineageFile" or die "Cannot open $newLineageFile for writing: $OS_ERROR";
my @tx;

foreach my $x (sort values %finaltx)
{
  push @tx, $x;
}
my @uniqTx = uniq @tx;

foreach my $values (@uniqTx)
{
  if ($lineage{$values})
  {
    print OUT "$values\t$lineage{$values}\n";
  }
  else
  {
    print "No lineage for $values\n";
  }
}

close OUT;
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

sub readArray{

  my $file = shift;
  my @ary;

  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  foreach (<IN>)
  {
    chomp;
    push @ary, $_;
  }
  close IN;

  return @ary;
}



exit;