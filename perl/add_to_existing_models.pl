#!/usr/bin/env perl

=head1 NAME

  add_to_existing_models.pl

=head1 DESCRIPTION

  Given the PECAN classification of a set of reference sequences, if the
  classification is at the species level, place the sequence ID and
  classified taxonomy into to_add.tx, and the sequence ID into to_add.seqID. 

  These files can then be used to replace particular species annotations in 
  the full dataset .tx file.

=head1 SYNOPSIS


=head1 OPTIONS

=over

=item B<--PECAN-file, -p>
  Input PECAN classification file. 

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

  add_to_existing_models.pl -p MC_order7_results.txt

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
  "PECAN-output|p=s"      => \my $pecanFile,
  "quiet"                 => \my $quiet,
  "verbose|v"             => \my $verbose,
  "debug"                 => \my $debug,
  "dry-run"               => \my $dryRun,
  "help|h!"               => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);

if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit 0;
}

####################################################################
##                               MAIN
####################################################################

my %pecan = read2colTbl($pecanFile);

my $txFile = "to_add.tx";
my $seqIDFile = "to_add.seqIDs";

open TX, ">$txFile" or die "Cannot open $txFile for writing: $OS_ERROR\n";
open SEQID, ">$seqIDFile" or die "Cannot open $seqIDFile for writing: $OS_ERROR\n";

my $count=0;
my $out=0;
foreach my $x (keys %pecan)
{
  $count++;
  my $sp = $pecan{$x};
  if ($sp =~ /^g_/ || $sp =~ /^f_/ || $sp =~ /^o_/ || $sp =~ /^c_/ || $sp =~ /^p_/ || $sp =~ /^d_/ )
  {
    next;
  }
  else
  {
    $out++;
    print TX "$x\t$pecan{$x}\n";
    print SEQID "$x\n";
  }
}

close TX;
close SEQID;

print "---Of $count sequences, $out sequences classified to the species level\n";

####################################################################
##                               SUBS
####################################################################

sub read2colTbl{

  my $file = shift;

  if ( ! -f $file )
  {
    warn "\n\nERROR in read2colTbl(): $file does not exist\n\n\n";
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

exit;