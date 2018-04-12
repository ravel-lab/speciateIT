#!/usr/bin/env perl

=head1 NAME

  score_classification.pl

=head1 DESCRIPTION

  Compares the expected versus the observed taxonomic classifications for a set of sequences 
  from CV runs considering the entire lineage of the expected taxonomy as follows:

    - Each species has 7 levels of taxonomy, so correct classification to each of those levels
      equals a point of 1.00/8 = 0.125 x the level reached. 
    - Correct classification to each level receives a score as follows: 
        Domain (d_):    1 x 0 = 0
        Phylum (p_):    2 x 0.125 = 0.25
        Class (c_):     3 x 0.125 = 0.325
        Order (o_):     4 x 0.125 = 0.5
        Family (f_):    5 x 0.125 = 0.625
        Genus (g_):     6 x 0.125 = 0.75
        Subgenus (sg_): 7 x 0.125 = 0.875
        Species (s_):   8 x 0.125 = 1

        Where 1 is the best score, and 0.125 is the worst score. 

    - Use the observed classification (following classify with error thresholds), and determine
      the accuracy of that classification.
    - The final score is printed out. A perfect score would equal nSeqs * 1 = nSeqs.
    - The per sequence scores are printed to the file: clScore.txt in the out directory provided.

=head1 SYNOPSIS

  score_classification.pl -i <input dir> -o <output dir>

  **Currently, looks for testSeqID.lineage and MC_order7_results.txt 
    as input from the in directory.

=head1 OPTIONS

=over

=item B<--input-dir, -p>
  PECAN output classification file.

=item B<--output-dir, -l>
  Lineage file with seqID => lineage, as produced by get_lineage.pl.

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

  score_classification.pl -i ~/Desktop/PECAN/FFT_NS_2/test1/cvDir/ -o ~/Desktop/PECAN/FFT_NS_2/test1/cvDir/

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

GetOptions(
  "lineage-file|l=s"    => \my $lineageFile,
  "PECAN-file|p=s"      => \my $pecanFile,
  "cluster-Tbl|c=s"     => \my $clusterTbl,
  "output-dir|o=s"      => \my $outDir,
  "verbose|v"           => \my $verbose,
  "debug"               => \my $debug,
  "dry-run"             => \my $dryRun,
  "help|h!"             => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);

my $testLineageFile;
my $clTxFile;

if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit 0;
}

if (!$lineageFile || !$pecanFile) 
{
  print "\nPlease provide a *_sppSeqID.lineage file (-l) and the PECAN classification file (-p)";
  exit 0;
}

####################################################################
##                               MAIN
####################################################################

my %testLineage = read2colTbl($lineageFile);
my %clTx = readTbl($pecanFile);
my @clIDs = keys %clTx;
my $test = $clIDs[1];
my %clScore;

if ($clusterTbl && -e $clusterTbl)
{
  my %cluster = read2colTbl($clusterTbl);
  %clScore = score_classification_w_cluster( \@clIDs, \%testLineage, \%clTx, \%cluster);
}
else
{
  %clScore = score_classification( \@clIDs, \%testLineage, \%clTx);
}

my $clScore_sum = sum values %clScore;
my $clScore_perc =  100 * $clScore_sum /  scalar @clIDs; 
my $rounded = sprintf "%.1f", $clScore_perc;
print "Mean overall classification score: " . $rounded;

my $nCorrect = 0; 
foreach my $s (values %clScore)
{
  if ($s == "1")
  {
    $nCorrect++;
  }
}
my $clScore_spp = $nCorrect / scalar @clIDs;
$rounded = sprintf "%.3f", $clScore_spp;
print "\nSpecies-level classification score: " . $rounded * 100;
my $clScoreFile = "clScore.txt";
my $clCompFile= "cl_cmp.txt";
if ($outDir)
{
  $clScoreFile = $outDir . "/clScore.txt";
  $clCompFile = $outDir . "/cl_cmp.txt";
}

open OUT, ">$clScoreFile" or die "Cannot open $clScoreFile for appending: $OS_ERROR";
open OUT1, ">$clCompFile" or die "Cannot open $clCompFile for appending: $OS_ERROR";
print OUT1 "seqID\tclassification\tcorrect_lineage_g\tcorrect_lineage_f\tcorrect_lineage_o\tcorrect_lineage_c\tcorrect_lineage_p\tcorrect_lineage_d\tcorrect_lineage_sp\n";
foreach my $key (keys %clScore)
{
  print OUT1 $key . "\t". $clTx{$key} ."\t". $testLineage{$key} . "\n";
  print OUT $key . "\t" . $clScore{$key} . "\n";

}
close OUT;
close OUT1;

print "\n---Per-sequence classification scores printed to $clScoreFile";
print "\n---Per-sequence classification comparison to reference printed to $clCompFile\n";

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
    if ($debug)
    {
    #print "\nMaking a hash table with $id connected to $tbl{$id}\n"
    }
  }
  close IN;

  return %tbl;
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

    ## t -> classified taxa
    my $t = $clTx{$s};
    if (!$t)
    {
      if ($debug)
      {
        print "\nWARNING: $s not found in $clTxFile\n";
      }
    } 
    else
    {
      ##print "The value of $s is $clTx{$s}.\n";
    }
    ## a -> reference lineage
    my @a = split("\t", $testLineage{$s});
    if (!@a)
    {
      if ($debug)
      {
        print "\nWARNING: $s not found in refLineage hash table\n;"
      }
      next;
    }
    ## r -> reference species
    my $r = $a[6];
    ## if the classified taxa is found in the reference lineage
    my $value = 0;
    if ( grep /$t/, @a) 
    {   
      if ($clTx{$s} =~ /^d_/) 
        {
          $clScore{$s} = 0;
        }
        elsif ($clTx{$s} =~ /^p_/) 
        {
          $clScore{$s} = 0.17;
        }
        elsif ($clTx{$s} =~ /^c_/) 
        {
          $clScore{$s} = 0.33;
        }
        elsif ($clTx{$s} =~ /^o_/) 
        {
          $clScore{$s} = 0.5;
        }
        elsif ($clTx{$s} =~ /^f_/) 
        {
          $clScore{$s} = 0.66;
        }
        elsif ($clTx{$s} =~ /^g_/) 
        {
          $clScore{$s} = 0.83;
        }          
        else 
        {
          $clScore{$s} = 1;
        }
    } 
    else ## And of course, these have to be correct with the lineage. If they are not correct, they get a zero.
    {
      $clScore{$s} = 0;
    } 
  }
print "Scores for ". scalar(keys %clScore) . " completed\n"; 
return %clScore;
}
sub score_classification_w_cluster{ 
  my $rclIDs = shift;
  my $rtestLineage = shift;
  my $rclTx = shift;
  my $rcluster = shift;

  my %testLineage = %$rtestLineage;
  my %clTx = %$rclTx;
  my %cluster = %$rcluster;
  my %clScore;

  my @clIDs = @$rclIDs;
  foreach my $s (@clIDs) ## For each classified sequence
  {
    chomp $s;
    ##print "\nHere is the lineage for $s: ". $testLineage{$s} . "\n";

    ## t is equal to the classified taxa
    my $class = $clTx{$s};
    my $refLin = $testLineage{$s};
    my $refsp;
    if (!$refLin)
    {
      if ($debug)
      {
        print "\nWARNING: $s not found in testLineage hash table\n;"
      }
      next;
    }
    else
    {
      my @lineage = split(/\t/, $refLin);
      #print Dumper \@lineage;
      $refsp = $lineage[6]; 
    }

    my $g;
    my $sp;
    #print "classification is $class\n";
    
    ($g, $sp) = split(/_/, $class);

    ## a is equal to the reference lineage
    
    #print "\nThe value of $s is $testLineage{$s}\n";
    ## if the classified taxa is found in the reference lineage

    my $value = 0;
    if ( grep /.*$class.*/, $refLin) 
    {   
      $value++;
      if ( $value > 0) ##If the classifed taxonomy is found in the array
      {
        if ($clTx{$s} =~ /^d_/) 
        {
          $clScore{$s} = 0;
        }
        elsif ($clTx{$s} =~ /^p_/) 
        {
          $clScore{$s} = 0.17;
        }
        elsif ($clTx{$s} =~ /^c_/) 
        {
          $clScore{$s} = 0.33;
        }
        elsif ($clTx{$s} =~ /^o_/) 
        {
          $clScore{$s} = 0.5;
        }
        elsif ($clTx{$s} =~ /^f_/) 
        {
          $clScore{$s} = 0.66;
        }
        elsif ($clTx{$s} =~ /^g_/) 
        {
          $clScore{$s} = 0.83;
        }          
        else 
        {
          $clScore{$s} = 1;
        }
      } 
      #print "\n1 The score for $s is $clScore{$s}\n";
    } 
    elsif ($class =~ /Cluster/)
    {
      $class = $cluster{$class};
      if ($class =~ /$refsp/ )
      {
        $clScore{$s} = 1;
      }
      #print "\n2 The score for $s is $clScore{$s}\n";
    }
    elsif ($refsp =~ /$class/)
    {
      $clScore{$s} = 1;
      #print "\n3 The score for $s is $clScore{$s}\n";
    }
    elsif ($class =~ /$refsp/)
    {
      $clScore{$s} = 1;
      #print "\n4 The score for $s is $clScore{$s}\n";
    }
    else ## And of course, these have to be correct with the lineage. If they are not correct, they get a zero.
    {
      $clScore{$s} = 0;
      #print "\n5 The score for $s is $clScore{$s}\n";
      if ($debug)
      {
        print "\n$class not found in $refLin\n";
      }
    }
    if ($debug)
    {
      print "$s is $refsp and was classified as $class and given score of $clScore{$s}\n"; 
    }
  }
print "Scores for ". scalar(keys %clScore) . " completed\n"; 
return %clScore;
}
exit 0;