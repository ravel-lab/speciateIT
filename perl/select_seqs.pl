#!/usr/bin/env perl

=head1 NAME

  select_seqs.pl

=head1  DESCRIPTION

  Select sequences from a single fasta file given a list of sequence
  IDs. Selected sequences are place in the new fasta file in the order they are
  listed.

=head1  INPUT

  - file with readIds of sequences to be selected
  - fasta file to select from

=head1  OUTPUT

  - fasta file of selected sequences

=head1 SYNOPSIS

    select_seqs.pl [-s|-e] <list of seqIDs file> -i <fastaFile> -o <outFile>

=head1 OPTIONS

=over

=item B<-h|--help>

    Print help message and exit successfully.

=item B<--fasta-file, -i>
  Fasta file

=item B<--sels-file, -s>
  File with a list of sequence IDs to be selected from the input fasta file

=item B<--excl-file, -e>
  File with a list of sequence IDs to be excluded from the input fasta file

=item B<--output-file, -o>
  Output file.

=item B<--quiet>
  Do not print progress messages.

=back

=head1  EXAMPLES

    select_seqs.pl -s Unclassified_nr.1k.sample.align.cmd.3d.cl2 -i Unclassified_nr.1k.sample.align.fa  -o Unclassified_nr.1k.sample.align.cmd.3d.cl2.fa

=cut

## Copyright (C) 2017 Pawel Gajer pgajer@gmail.com
##
## Permission to use, copy, modify, and distribute this software and its
## documentation with or without modifications and for any purpose and
## without fee is hereby granted, provided that any copyright notices
## appear in all copies and that both those copyright notices and this
## permission notice appear in supporting documentation, and that the
## names of the contributors or copyright holders not be used in
## advertising or publicity pertaining to distribution of the software
## without specific prior permission.
##
## THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
## WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
## WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
## CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
## OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
## OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
## OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
## OR PERFORMANCE OF THIS SOFTWARE.
##


use strict;
use warnings;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

GetOptions(
  "input-file|i=s"  => \my $inFile,
  "sels-file|s=s"   => \my $selsFile,
  "excl-file|e=s"   => \my $exclFile,
  "output-file|o=s" => \my $outFile,
  "quiet"           => \my $quiet,
  "help|h!"         => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);


if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if (!$inFile)
{
  print "\n\nERROR: Missing input file\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif (!$outFile)
{
   print "\n\nERROR: Missing output file\n\n\n";
   pod2usage(verbose => 2,exitstatus => 0);
   exit 1;
}

if ( ! -e $inFile )
{
  print "\n\nERROR: $inFile does not exist\n\n\n";
  exit 1;
}

if ( -l $inFile )
{
  $inFile = readlink($inFile);
}

if ( $selsFile && -l $selsFile )
{
   $selsFile = readlink($selsFile);
}

if ( $exclFile && -l $exclFile )
{
  $exclFile = readlink($exclFile);
}

####################################################################
##                               MAIN
####################################################################

my $count = 0;
my $selCount = 0;

if ( defined $selsFile )
{
  my %selRead = read_array($selsFile);
  my $nSelSeqs = keys %selRead;

  open (FASTA, "<$inFile") or die "Cannot open $inFile for reading: $OS_ERROR\n";
  $/ = ">";
  my $junkFirstOne = <FASTA>;

  open OUT, ">$outFile" or die "Cannot open $outFile for reading: $OS_ERROR\n";
  while (<FASTA>)
  {
    chomp;
    my ($def,@seqlines) = split /\n/, $_;
    my $seq = join '', @seqlines;
    my ($id) = split /\s+/, $def;
    if ( exists $selRead{$id} && $selRead{$id}==1 )
    {
      print OUT ">$id\n$seq\n";
      $selCount++;
      $selRead{$id}++;
    }
    last if $selCount==$nSelSeqs;
    $count++;
  }
  close OUT;
  close FASTA;

  #print "nSelSeqs: $nSelSeqs\n";
  my @missed;
  for ( keys %selRead )
  {
    if ( $selRead{$_}==1 )
    {
      push @missed, $_;
    }
  }
  if ( @missed )
  {
    print "\n\nSequence IDs not found in the fasta file\n";
    print_array( \@missed );
    print "\n";
  }
}
else
{
  my %exclRead = read_array($exclFile);

  open (FASTA, "<$inFile") or die "Cannot open $inFile for reading: $OS_ERROR\n";
  $/ = ">";
  my $junkFirstOne = <FASTA>;
  open OUT, ">$outFile" or die "Cannot open $outFile for reading: $OS_ERROR\n";
  while (<FASTA>)
  {
    chomp;
    my ($def,@seqlines) = split /\n/, $_;
    my $seq = join '', @seqlines;
    my ($id) = split /\s+/, $def;
    if ( !exists $exclRead{$id} )
    {
      print OUT ">$id\n$seq\n";
      $selCount++;
    }
    $count++;
  }
  close OUT;
  close FASTA;
  $/ = "\n";
}

if (!$quiet)
{
  print "\r                                                           ";
  print "\n\tNumber of sequences in the input fasta file:  $count\n";
  print "  \tNumber of sequences in the output fasta file: $selCount\n\n";

  print "\rSelected sequences written to $outFile\n";
}

####################################################################
##                               SUBS
####################################################################

# print array to stdout
sub print_array
{
  my ($a, $header) = @_;
  print "$header\n" if $header;
  map {print "$_\n"} @{$a};
}

# read table with one column
sub read_array
{
  my $file = shift;
  my %rows;

  #local $/ = '\n';
  open (INPUT, "<", $file) or die "Cannot open $file for reading: $OS_ERROR\n";
  foreach (<INPUT>)
  {
    chomp;
    my ($id) = split /\s+/, $_;
    #print "id: $id\n"; exit 1;
    $rows{$id}=1;
  }
  close INPUT;

  return %rows;
}

# print elements of a hash table
sub printTbl{

  my $rTbl = shift;
  map {print "|$_|\t|" . $rTbl->{$_} . "|\n"} keys %$rTbl;
}

# print n keys of a hash table
sub printKeys{

  my ($rTbl, $n) = @_;
  my $i = 0;
  map {$i++; print "$_\n" if $i < $n} keys %$rTbl;
}

exit 0;
