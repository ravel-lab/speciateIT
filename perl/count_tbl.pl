#!/usr/bin/perl

=head1 NAME

  count_tbl.pl

=head1 DESCRIPTION

  given a taxon file with two columns <read ID> <OTU>
  generate a corresponding sampleID x OTU contingency table
  assuming that each read ID is of the form

  subjID_visit_readIndex
  OR
  subjID.visit_readIndex

=head1 SYNOPSIS

  count_tbl.pl -i <taxon file> -o <output file> [Options]

=head1 OPTIONS


=over

=item B<--taxon-file, -i>
  Taxon file.

=item B<--output-file, -o>
  Output file.

=item B<--dry-run>
  Print commands to be executed, but do not execute them.

=item B<-h|--help>
  Print help message and exit successfully.

=back


=head1 EXAMPLE

  cd /usr/local/projects/pgajer/projects/microbicide/rdpTxDir/Lactobacillus

  count_tbl.pl -i Lactobacillus.taxonomy -o Lactobacillus.spp.count.tbl.txt

=cut

use strict;
use warnings;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use File::Basename;
use List::Util qw(min);

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

GetOptions(
  "taxon-file|i=s"    => \my $txFile,
  "output-file|o=s"   => \my $outFile,
  "dry-run"           => \my $dryRun,
  "help|h!"           => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);


if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}

if (!$txFile)
{
  print "ERROR: Missing input file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}
elsif (!$outFile)
{
  print "ERROR: Missing output file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}

####################################################################
##                               MAIN
####################################################################


my $startRun = time();
my $endRun = time();
my $runTime = $endRun - $startRun;
my $timeStr;
my $timeMin = int($runTime / 60);
my $timeSec = $runTime % 60;

print "--- Parsing $txFile";
my %txTbl = read2colTbl($txFile); # <seq Id> => <taxonomy>

print "\r--- Creating sample x OTU count table";
my %sampleTxTbl;
my %sampleTbl;
my %otuTbl;
my %rankTbl;# rank table
my $rcount = 0;
my $counter = 1;
my $nTx = keys %txTbl;
##my $newSample = 0;
for my $readId (keys %txTbl)
{
   if ($counter % 500 == 0)
   {
      $endRun = time();
      $runTime = $endRun - $startRun;
      if ( $runTime > 60 )
      {
         $timeMin = int($runTime / 60);
         $timeSec = sprintf("%02d", $runTime % 60);
         $timeStr = "$timeMin:$timeSec";
      }
      else
      {
         $runTime = sprintf("%02d", $runTime);
         $timeStr = "$timeMin:$runTime";
      }
      my $perc = sprintf("%.1f%%", 100 * $counter / $nTx);
      print "\r$timeStr [$perc]";
   }
   $counter++;

   # my ($id) = split ":", $readId;
   # my ($sampleID, $index) = ( $id =~ /(.+)_(\d+$)/ );

   # if ( !exists $sampleTbl{$sampleID} )
   # {
   #   $newSample++;
   #   print "\n\nsample $newSample\treadID: $readId\n";
   #   print "sampleID: $sampleID\tindex: $index\n";
   # }

   my $otu = $txTbl{$readId};
   print "ERROR: otu undefined for $readId\n" if !$otu;

   # print "$readId => $otu\n";

   $otuTbl{$otu} = 0;
   $rankTbl{$otu}++;

   my $sampleID = $readId;

   if ($sampleID)
   {
      $sampleTbl{$sampleID}++;
      $sampleTxTbl{$sampleID}{$otu}++;
   }
   else
   {
      print "\n\n==== ERROR: in " . __FILE__ . " at line " . __LINE__ . ": $OS_ERROR\n";
      print "\nreadId: $readId\n";
      ##print "id: $id\n\n";
   }

   $rcount++;
}

## compute column sums for sorting OTUs
print "\r--- Sorting OTUs by column sums";
my %colSums;
my @tIds = keys %otuTbl;
my $ntIds = @tIds;

my $i = 0;
foreach my $tid (@tIds)
{
   print "\r$i\t$tid";
   $i++;
   my $sum = 0;
   foreach my $sid (keys %sampleTxTbl)
   {
      $sum += $sampleTxTbl{$sid}{$tid} if exists $sampleTxTbl{$sid}{$tid};
   }
   $colSums{$tid} = $sum;
}
print "\r                                                                                                                         \n";

my $nTids = @tIds;
my @otus = sort { $colSums{$b} <=> $colSums{$a} } @tIds;


## computing max string length of the first $n elements of @otus

my $n = min(scalar(@tIds), 20);
my $n1 = $n - 1;
print "\n\tFrequencies of the $n most abundant phylotypes\n\n";

my %maximums;
$maximums{"riboID"} = 0;
$maximums{"count"} = 0;
$maximums{"perc"} = 0;

my @otusR = @otus[0..$n1];

# calculate the maximum length of the values in each column
my @fTbl;
foreach my $rID (@otusR)
{
   my %row;
   $row{"riboID"} = $rID;
   $row{"count"}  = commify($rankTbl{$rID});
   $row{"perc"}   = sprintf("%.2f", (100*$rankTbl{$rID}/$rcount));
   push @fTbl, \%row;
   foreach my $key (keys %maximums)
   {
      my $col_length = length($row{$key});
      $maximums{$key} = $col_length if ($col_length > $maximums{$key});
   }
}

## adding extra space
my $dx = 5;
foreach my $key (keys %maximums)
{
   $maximums{$key} += $dx;
}

## format string
my $row_format = "%" . $maximums{"riboID"} . "s";
$row_format   .=  "  %" . $maximums{"count"}  . "s";
$row_format   .=  "  %" . $maximums{"perc"}  . "s\n";

printf($row_format, "Phylotype", "Frequency", "Percentage");
foreach my $row (@fTbl) {
   printf($row_format, $row->{"riboID"}, $row->{"count"}, $row->{"perc"});
}
print "\n\n";


## print "\r                                     \nRank\tFrequencies\tPercentages\n";
## map {print "$_\t" . $rankTbl{$_} . "\t" . (100*$rankTbl{$_}/$rcount) . "\n"} keys %rankTbl;
## map {print "$_\t" . sprintf("%10s", commify($rankTbl{$_})) . "\t" . sprintf("%10.6f\n", (100*$rankTbl{$_}/$rcount)) } @otus[0..$n1];

my $nOTUs = @otus;
my $nSamp = keys %sampleTbl;

print "\n\tTotal number of seq's: " . sprintf("%10s", commify($rcount)) . "\n";
print   "\tNumber of phylotypes:  " . sprintf("%10s", commify($nOTUs)) . "\n";
print   "\tNumber of samples:     " . sprintf("%10s", commify($nSamp)) . "\n\n";

my $txRkFile = dirname($outFile) . "/taxonRank.txt";
open OUT, ">$txRkFile" or die "Cannot open $txRkFile for writing: $OS_ERROR\n";
map {print OUT "$_\t" . $rankTbl{$_} . "\t" . (100*$rankTbl{$_}/$rcount) . "\n"} @otus;
close OUT;

print "\n\tTaxon frequencies written to $txRkFile\n";

my %otuIdx;
$i = 0;
foreach my $g (@otus)
{
   $otuIdx{$g} = $i;
   $i++;
}

## print "--- Printing sample x OTU count table to $outFile\n";
open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR";
# header
print OUT "sampleID";
foreach my $otu (@otus)
{
   print OUT "\t$otu";
}
print OUT "\n";

$counter = 1;
my $nStbl = keys %sampleTbl;
foreach my $sid (sort keys %sampleTbl)
{
   if ($counter % 500 == 0)
   {
      $endRun = time();
      $runTime = $endRun - $startRun;
      if ( $runTime > 60 )
      {
         $timeMin = int($runTime / 60);
         $timeSec = sprintf("%02d", $runTime % 60);
         $timeStr = "$timeMin:$timeSec";
      }
      else
      {
         $timeStr = "$timeMin:$runTime";
      }
      my $perc = sprintf("%.1f%%", 100 * $counter / $nStbl);
      print "\r$timeStr [$perc]";
   }
   $counter++;

   my @counts = (0) x @otus;

   foreach my $otu (keys %{$sampleTxTbl{$sid}})
   {
      $counts[$otuIdx{$otu}] = $sampleTxTbl{$sid}->{$otu};
   }

   print OUT $sid;
   foreach my $count (@counts)
   {
      print OUT "\t$count";
   }
   print OUT "\n";
}
close OUT;

print "\r\tSample x Phylotype count table written to $outFile\n\n";

## report timing
if (0)
{
   $endRun = time();
   $runTime = $endRun - $startRun;
   if ( $runTime > 60 )
   {
      $timeMin = int($runTime / 60);
      $timeSec = sprintf("%02d", $runTime % 60);
      print "\rCompleted in $timeMin:$timeSec\n"
   }
   else
   {
      print "\rCompleted in $runTime seconds\n"
   }
}



####################################################################
##                               SUBS
####################################################################

## put commas in numbers for better readability
## lifted from
## http://www.perlmonks.org/?node_id=2145
sub commify {
   local $_  = shift;
   s{(?<!\d|\.)(\d{4,})}
   {my $n = $1;
     $n=~s/(?<=.)(?=(?:.{3})+$)/,/g;
     $n;
    }eg;
   return $_;
}

# read two column table; create a table that assigns
# elements of the first column to the second column
sub read2colTbl{

   my $file = shift;

   ## my $fileSizeL = -s $file; # file size in bytes
   ## my %vals;
   my %tbl;
   my $counter = 1;
   open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
   foreach (<IN>)
   {
      if ($counter % 500 == 0)
      {
         $endRun = time();
         $runTime = $endRun - $startRun;
         if ( $runTime > 60 )
         {
            $timeMin = int($runTime / 60);
            $timeSec = sprintf("%02d", $runTime % 60);
            $timeStr = "$timeMin:$timeSec";
         }
         else
         {
            $runTime = sprintf("%02d", $runTime);
            $timeStr = "$timeMin:$runTime";
         }
         ##my $perc = sprintf("%.1f%%", 100 * (tell IN) / $fileSizeL);
         ##print "\r$timeStr [$perc]";
         print "\r$timeStr";
      }
      $counter++;

      chomp;
      my ($id, $t) = split /\s+/,$_;
      $tbl{$id} = $t;

      # if ( !exists $tbl{$id} )
      # {
      #   $tbl{$id} = $t;
      # }
      # else
      # {
      #   print "rec: $_\n";
      #   print "id: $id\ttx: $t\ttbl{$id}: " . $tbl{$id} . "\n";
      #   exit;
      # }
      ## $vals{$t} = 0;
   }
   close IN;

   return %tbl;
}

# extract unique elements from an array
sub unique{

   my $a = shift;
   my %saw;
   my @out = grep(!$saw{$_}++, @{$a});

   return @out;
}

# print array to stdout
sub printArray{

   my ($a, $header) = @_;

   print "$header\n" if $header;
   map {print "$_\n"} @{$a};
}

# print elements of a hash table
sub printTbl{

   my $rTbl = shift;
   map {print "$_\t" . $rTbl->{$_} . "\n"} keys %$rTbl;
}

# format OUT name
# sub formatOTUname{

#   my $otu = shift;

#   if ($otu =~ /^c/)
#   {
#     $otu =~ s/^c/L.otu/;
#     $otu =~ s/\.\d+/$idx/;
#     $idx++;
#   }
#   else
#   {
#     $otu =~ s/\.\d+//;
#   }

#   return $otu;
# }

exit;
