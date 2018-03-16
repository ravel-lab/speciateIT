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

=item B<--input-dir, -i>
  Input directory containing all.fa, spp.tx, and spp.lineage files. 

=item B<--output-dir, -o>
  Output directory.

=item B<--input-dir, -c>
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
  "input-dir|i=s"       => \my $inDir,
  "output-dir|o=s"      => \my $outDir,
  "lineage-file|l=s"    => \my $l,
  "MC-file|c=s"         => \my $c,
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

if ( !$inDir && !$l && !$c)
{
  warn "\n\n\tERROR: No input directory or files provided.";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 0;
}
elsif ( !$inDir && !$l )
{
  warn "\n\n\tERROR: No input directory or lineage file provided.";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 0;
}
elsif ( !$inDir && !$c )
{
  warn "\n\n\tERROR: No input directory or classification file provided.";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 0;
}
elsif ( !$inDir )
{
  print "\n\n\tNo input directory provided.";
  print "\n\tUsing $l and $c files.";
  $testLineageFile = $l;
  $clTxFile = $c;
}
elsif (!$l && !$c) 
{
  $testLineageFile = $inDir . "/sppSeqID.lineage";
  $clTxFile = $inDir . "/MC_order7_results.txt";
}

if ( $inDir && !$outDir )
{
  print "\n\n\tERROR: Missing out directory. Using input directory for output if provided.";
  print "\n\n";
  $outDir = $inDir;
}
elsif ( !$inDir && !$outDir )
{
  $outDir = getcwd;
}



####################################################################
##                               MAIN
####################################################################

my %testLineage = read2colTbl($testLineageFile);
my %clTx = readTbl($clTxFile);
my @clIDs = keys %clTx;
my $test = $clIDs[1];

my %clScore = score_classification( \@clIDs, \%testLineage, \%clTx);
my $clScore_sum = sum values %clScore;
my $clScore_perc =  100 * $clScore_sum /  scalar keys %clScore; 
print "\n\nThe classification score for ". scalar @clIDs . " sequences is: " . $clScore_perc . "\n";
print "Here is a value of the clScore hash table: " . $clScore{$test} ."\n\n";
my $clScoreFile = "$outDir/clScore.txt";
my $clCompFile= "$outDir/cl_cmp.txt";
  open OUT, ">$clScoreFile" or die "Cannot open $clScoreFile for appending: $OS_ERROR";
  open OUT1, ">$clCompFile" or die "Cannot open $clCompFile for appending: $OS_ERROR";
  print OUT1 "seqID\tclassified_lineage\tcorrect_lineage_sg\tcorrect_lineage_g\tcorrect_lineage_f\tcorrect_lineage_o\tcorrect_lineage_c\tcorrect_lineage_p\tcorrect_lineage_d\tcorrect_lineage_sp\n";
  foreach my $key (keys %clScore)
  {
    print OUT1 $key . "\t". $clTx{$key} ."\t". $testLineage{$key} . "\n";
    print OUT $key . "\t" . $clScore{$key} . "\n";

  }
  close OUT;
  close OUT1;

@clIDs = keys %clScore;

print "\nClassification scores printed to $clScoreFile\n";
print "Classification comparison to reference printed to $clCompFile\n\n";
print "\nThe classification score for ". scalar @clIDs . " sequences is: " . $clScore_perc . "\n";


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

    ## t is equal to the classified taxa
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
return %clScore;
}
exit 0;