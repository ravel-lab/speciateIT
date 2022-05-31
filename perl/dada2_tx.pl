#!/usr/bin/env perl

=head1 NAME

  dada2_tx.pl

=head1 DESCRIPTION

  Collapsing dada2 and PECAN taxonomy of data2 tables into a phylotype count
  table. Creating also a taxonomy row in the data2 table.

=head1 SYNOPSIS

  dada2_tx.pl -i <dada2 file> -t <taxonomy of ASV's> -o <output dir> [Options]

=head1 OPTIONS

=over

=item B<--dada-file, -i>
  CSV dada2 file.

=item B<--tx-file, -t>
  Tab delimited taxonomy file.

=item B<--output-dir, -o>
  Output directory.

=item B<--verbose, -v>
  Prints content of some output files.

=item B<-h|--help>
  Print help message and exit successfully.

=back


=head1 EXAMPLE

  dada2_tx.pl -i all_runs_dada2_abundance_table_all_controls.csv -t MC_order7_results_all_controls.txt -o ctls_dir

=cut

use strict;
use warnings;
use diagnostics;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use List::Util qw( sum );

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

GetOptions(
   "dada-file|i=s"  => \my $dadaFile,
   "tx-file|t=s"    => \my $txFile,
   "out-dir|o=s"    => \my $outDir,
   "verbose|v"      => \my $verbose,
   "help|h!"        => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);


if ( $help )
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if ( !$dadaFile )
{
  print "\n\n\tERROR: Missing dada2 file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif ( !$txFile )
{
  print "\n\n\tERROR: Missing taxonomy file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif ( !$outDir )
{
  print "\n\n\tERROR: Missing output dir\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

####################################################################
##                               MAIN
####################################################################

print "--- Parsing taxonomy file\n";
my %txTbl = parse_tx( $txFile );

print "--- Parsing dada2 file\n";
my ($rtbl, $rttbl, $rrowIds, $rheader) = parse_dada( $dadaFile );

my $cmd = "rm -rf $outDir; mkdir $outDir";
system($cmd) == 0 or die "system($cmd) failed:$?";

print "--- Creating dada2 file with taxonomy\n";
my $outFile = $outDir . "/dada2_tx.csv";
my %dadaTbl = %{$rtbl};
open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
## ASV's
map{ print OUT ",$_" } @{$rheader};
print OUT "\n";
## their taxonomy
map{ print OUT "," . $txTbl{$_} } @{$rheader};
print OUT "\n";
## counts
for my $id ( keys %dadaTbl )
{
   print OUT "$id";
   for ( @{$dadaTbl{$id}} )
   {
      print OUT ",$_";
   }
   print OUT "\n";
}
close OUT;


print "--- Creating phylotype count table\n";

## inversing taxonomy table
my %invTx;
for ( keys %txTbl )
{
   push @{ $invTx{$txTbl{$_}} }, $_;
}

my %tDadaTbl = %{$rttbl};
my $nSamples = keys %dadaTbl;
my %phTrTbl;
for my $tx ( keys %invTx )
{
   my @counts = (0) x $nSamples;
   for my $avs ( @{$invTx{$tx}} )
   {
      ##add_to( $ttbl{$avs}, \@counts );
      @counts = map{ $tDadaTbl{$avs}->[$_] + $counts[$_] } 0..$#counts;
   }
   $phTrTbl{$tx} = \@counts;
}

# transposing phTrTbl

my @txs = sort { sum( @{$phTrTbl{$b}} ) <=> sum( @{$phTrTbl{$a}} ) } keys %phTrTbl;

$outFile = $outDir . "/phylotype_read_count_tbl.csv";
open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
map{ print OUT ",$_" } @txs;
#my %phTbl;
my @sampleIDs = @{$rrowIds};
for my $i ( 0..$#sampleIDs )
{
   print OUT $sampleIDs[$i];
   for my $tx ( @txs )
   {
      ##push @{$phTbl{$sampleIDs[$i]}}, $phTrTbl{$tx}->[$i];
      print OUT "," . $phTrTbl{$tx}->[$i];
   }
   print OUT "\n";
}
close OUT;

print "Output written to $outDir\n";

####################################################################
##                               SUBS
####################################################################

# add an array to another array
# sub add_to
# {
#    my ($a1, $a2) = @_;
# }

# parse taxonomy table
sub parse_tx
{
  my $file = shift;

  if ( ! -e $file )
  {
    warn "\n\n\tERROR in parse_tx(): $file does not exist";
    print "\n\n";
    exit 1;
  }

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
  foreach (<IN>)
  {
    chomp;
    my ($id, $t) = split /\s+/,$_;
    $tbl{$id} = $t;
  }
  close IN;

  return %tbl;
}


# parse table with row ids (and optional header)
# output written to
# 1. tbl: <rowId> => <row of the table associated with rowId>
# 2. ttbl (transposed tbl): <colId> => <column of the table associated with colId>
# 3. rowIds: array of row ids, in the order as they appear in the input file
# 4. header array: if header present
#
sub parse_dada
{
   my $file = shift;

   if ( ! -e $file )
   {
      warn "\n\n\tERROR in parse_dada(): $file does not exist";
      print "\n\n";
      exit 1;
   }

   my %tbl;
   my @header;
   my @rowIds;

   open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";

   my $headerStr = <IN>;
   chomp $headerStr;
   $headerStr =~ tr/"//d;
   @header = split ",", $headerStr;
   shift @header; # get rid of <rowId> label

   foreach my $line (<IN>)
   {
      chomp $line;
      my @f = split ",", $line;
      my $id = shift @f;
      push @rowIds, $id;
      $tbl{ $id } = \@f;
   }
   close IN;

   my $ncol = @header;
   my $nrow = @rowIds;

   # transposed table
   my %ttbl;
   for my $i (0..$#header)
   {
      for ( keys %tbl )
      {
         push @{ $ttbl{$header[$i]} }, $tbl{$_}[$i];
      }
   }

  return (\%tbl, \%ttbl, \@rowIds, \@header);
}

# write hash table to a file
sub write_tbl
{
  my ($rTbl, $outFile) = @_;
  my %tbl = %{$rTbl};
  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR";
  map {print OUT $_ . "\t" . $tbl{$_} . "\n"} sort keys %tbl;
  close OUT;
}

# extract unique elements from an array
sub unique{

  my $a = shift;
  my %saw;
  my @out = grep(!$saw{$_}++, @{$a});

  return @out;
}

# common part of two arrays
sub comm
{
  my ($a1, $a2) = @_;

  my @c; # common array
  my %count;

  foreach my $e ( unique($a1), unique($a2) ){ $count{$e}++ }

  foreach my $e (keys %count)
  {
    push @c, $e if $count{$e} == 2;
  }

  return @c;
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

## Testing if two arrays are identical in a set-theoretic sense. That is that
## they have exactly the same set of elements.
sub set_equal
{
  my ($rA, $rB) = @_;

  my @a = @{$rA};
  my @b = @{$rB};
  my @c = comm(\@a, \@b);

  my $ret = 1;

  if (@c != @a || @c != @b)
  {
    warn "\n\n\tERROR: Elements of the two arrays do not match";
    print "\n\tNumber of elements in the first array: " . @a . "\n";
    print "\tNumber of elements in the second array: " . @b . "\n";
    print "\tNumber of common elements: " . @c . "\n";

    # writeArray(\@a, "a.txt");
    # writeArray(\@b, "b.txt");
    #print "\n\tNew taxon keys and fasta IDs written to a.txt and b.txt, respectively\n\n";

    if (@a > @b)
    {
      my @d = diff(\@a, \@b);
      print "\nElements a but not b:\n";
      for (@d)
      {
         print "\t$_\n";
      }
      print "\n\n";
    }

    if (@b > @a)
    {
      my @d = diff(\@b, \@a);
      print "\nElements in b that are not a:\n";
      for (@d)
      {
         print "\t$_\n";
      }
      print "\n\n";
    }

    $ret = 0;
  }

  return $ret;
}

exit 0;
