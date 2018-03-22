#!/usr/bin/env perl

=head1 NAME

  get_lineage.pl

=head1 DESCRIPTION

  Given a taxonomy file, and a lineage file, produce an output
  lineage file with seqID -> lineage
  
=head1 SYNOPSIS

  get_lineage.pl -t <taxonomic annotation file> -l <source lineage file>

=head1 OPTIONS

=over

=item B<--taxonomic-annotation file, -t>
  List of selected outgroup sequence IDs. 

=item B<--source-lineage-file, -l>
  The starting, source file containing lineage information for all sequences.
 
=item B<--output file, -o>
 Out file name.

=item B<--verbatim, -v>
  Prints content of some output files.

=item B<--debug>
  Prints system commands

=item B<-h|--help>
  Print help message and exit successfully.

=back

=head1 EXAMPLE


  
=cut

use strict;
use warnings;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Data::Dumper qw(Dumper);
use File::Basename;

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

GetOptions(
  "taxonomy-file|t=s"     => \my $txFile,
  "lineage-file|l=s"      => \my $lineageFile,
  "out-file|o=s"          => \my $outFile,
  "verbose|v"             => \my $verbose,
  "debug"                 => \my $debug,
  "dry-run"               => \my $dryRun,
  "help|h!"               => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);

if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}

if ( !$txFile )
{
  print "\n\nERROR: Missing taxonomic annotation file.\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}
elsif (!$lineageFile)
{
  print "\n\nERROR: Missing source lineage file.\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}
if (!$outFile)
{
    print "\nOutput file name not provided. Writing to new_sppSeqID.lineage.\n";
    $outFile = "new_sppSeqID.lineage";
}


####################################################################
##                               MAIN
####################################################################

## From a list of sequences get the full taxonomic lineage from the source file.
my $taxa;
my $lineage;

my %source = read2colTbl_and_add_column($lineageFile);
my %tx = readTbl($txFile);
my %done;

my $newLineage = $outFile;

open (OUT, ">$newLineage") or die "Cannot open $newLineage for reading: $OS_ERROR\n";
for my $seqID (keys %tx) 
{
  $taxa = $tx{$seqID};
  $lineage = $source{$taxa};
  if (!$lineage)
  {
    print "WARNING: $taxa does not exist in $lineageFile.\n"
  }
  elsif (!$taxa)
  {
    print "WARNING: $seqID does not exist in $txFile.\n"
  }
  else 
  {
    print OUT "$seqID\t$lineage\n";
  }
}
close OUT;

print "\nLineages for ". scalar(keys %tx) . " sequences from $txFile written to $newLineage.\n\n";

####################################################################
##                               SUBS
####################################################################


sub readTbl{

  my $file = shift;

  if ( ! -f $file )
  {
    print "\n\nERROR: $file does not exist\n\n\n";
    exit;
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

# read two column table
sub read2colTbl_and_add_column{

  my $file = shift;

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  foreach (<IN>)
  {
    chomp;
    next if $_ eq "";
    my ($id, $t) = split (/\s+/,$_, 2);
    $tbl{$id} = $t;
    $tbl{$id} = $tbl{$id} . "\t" . $id;
    if ($debug)
    {
    print "\nMaking a hash table with $id connected to $tbl{$id}\n"
    }
  }
  close IN;

  return %tbl;
}

exit;
