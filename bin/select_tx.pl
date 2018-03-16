#! /usr/bin/perl

=head1 NAME

  select_tx.pl

=head1 DESCRIPTION

  select sequence taxon assignment rows for selected seq IDs

=head1 SYNOPSIS

  select_tx.pl [-s|-e] <seq IDs> -i <input file> -o <output file> [Options]

=head1 OPTIONS


=over

=item B<--input-file, -i>
  Input file.

=item B<--seq-ids, -s>
  Sequence ids to be selected

=item B<--excl-file, -e>
  A file with a list of sequence IDs to be excluded from the input file

=item B<--output-file, -o>
  Output file.

=item B<--quiet>
  Do not print progress messages.

=item B<-h|--help>
  Print help message and exit successfully.

=back


=head1 EXAMPLE

  cd /Users/pgajer/devel/speciateIT/Lb_data

  select_tx.pl -s L.crispatus_trimmed_89_625_nr.ids -i crispatus.tx -o L.crispatus_trimmed_89_625_nr.tx

=cut

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
  "output-file|o=s" => \my $outFile,
  "seq-ids|s=s"     => \my $selsFile,
  "excl-file|e=s"   => \my $exclFile,
  "help|h!"         => \my $help,
  "quiet"           => \my $quiet,
  )
  or pod2usage(verbose => 0,exitstatus => 1);


if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if ( $help || !$inFile || !$outFile )
{
  pod2usage( { -exitval => 0, -verbose => 2 } );
  exit 1;
}

if ( !$selsFile && !$exclFile )
{
  print "Either -s or -e option has to be used.\n";
  pod2usage( { -exitval => 0, -verbose => 2 } );
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

my %tx = readTxTbl($inFile);

if ( $selsFile )
{
  my %seqIDs = readArray($selsFile);

  open OUT, ">$outFile" or die "Cannot open $outFile for reading: $OS_ERROR\n";
  for my $id (keys %tx)
  {
    if (exists $seqIDs{$id})
    {
      #print OUT $_ . "\t" . $tx{$_} . "\n";
      print OUT $id;
      map { print OUT "\t" . $_  }  @{$tx{$id}};
      print OUT "\n";
    }
  }
  close OUT;
}
elsif ( $exclFile )
{
  my %exclIDs = readArray($exclFile);

  open OUT, ">$outFile" or die "Cannot open $outFile for reading: $OS_ERROR\n";
  for my $id (keys %tx)
  {
    if (! exists $exclIDs{$id})
    {
      #print OUT $_ . "\t" . $tx{$_} . "\n";
      print OUT $id;
      map { print OUT "\t" . $_  }  @{$tx{$id}};
      print OUT "\n";
    }
  }
  close OUT;
}
else
{
  print "ERROR: no selection nor exclusion files\n\n";
  exit 1;
}

print "Output written to $outFile\n" if !$quiet;

####################################################################
##                               SUBS
####################################################################

# read table with one column
sub readArray{

  my ($file) = @_;
  my %rows;

  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  foreach (<IN>)
  {
    chomp;
    next if $_ eq "";
    my ($id) = split /\s+/,$_;
    $rows{$id}=0;
  }
  close IN;

  return %rows;
}

# read two column table; create a table that assigns
# elements of the first column to the second column
sub readTxTbl{

  my $file = shift;

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  foreach (<IN>)
  {
    chomp;
    my @t = split /\s+/,$_;
    my $id = shift @t;
    $tbl{$id} = \@t;
  }
  close IN;

  return %tbl;
}

exit 0;
