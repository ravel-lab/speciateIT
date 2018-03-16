my $rseqs_to_merge = shift;
my @seqs_to_merge = @$rseqs_to_merge;
my $rnewtx = shift;
my %newtx = %$rnewtx;
my $rlineage = shift;
my %lineage = %$rlineage;
my $rpecan = shift;
my %pecan = %$rpecan;

my %doneSeqs;

foreach my $seq (@seqs_to_merge)
{
  my @clSeqs;
  my $winner;
  my $newsp;
  my $refsp;
  my $clsp;
  my @refseqID;
  my @clseqID;
  my $nrefseq;
  my $nclseq;
  my $r;
  my $c;

  if (!exists $doneSeqs{$seq})
  {
 
    ## For each sequence, identify the species it was classified as and the species it came from
    ## compare the number of sequences for each of these species
    ## and which ever one has more sequences, change the lineage (ph - sg) to the one with more seq's
    ## and change the species name for all to contenate as biggersp_smallersp
    ## Therefore, we need the species it was classified as (pecanFile), the reference (%tx), the number
    ## of seq's in each of those species (from the ref %tx)

    $refsp = $newtx{$seq};
    @refseqID = grep { $tx{$_} eq $refsp } keys %newtx;
    print "There are " . scalar @refseqID . " sequences in $refsp\n";
    $nrefseq = scalar @refseqID;

    $clsp = $pecan{$seq};
    @clseqID = grep { $newtx{$_} eq $clsp } keys %newtx;
    print "There are " . scalar @clseqID . " sequences in $clsp\n";
    $nclseq = scalar @clseqID;

    $r = scalar @refseqID;
    $c = scalar @clseqID;

    if ( $r > $c)
    {
      $winner = 0;
    }
    elsif ( $r < $c)
    {
      $winner = 1;
    }
    elsif ( $r == $c)
    {
      $winner = 2;
      print "\n***WARNING: There are the same number of seq's in $refsp and $clsp\n";
    }
    ## These are all seqs in mc_output classified to this species.
    @clSeqs = grep { $pecan{$_} eq $clsp } keys %pecan;
      ## If the reference sp has more seq's, then we need to 
      ## 1. update the taxonomy (seqid -> tx) to be refsp_clsp
      ## 2. delete the clsp from the lineage hash
      ## 3. add a new key/value with the refsp_clsp to the ref_lineage
    if ($winner == 0)
    {
      $newsp = "$refsp"."_"."$clsp";
      chomp $newsp;
      if ($debug)
      {
        print "$clsp is now $newsp\n";
      }
      foreach my $clseq ( @clSeqs)
      {
        {
          $newtx{$clseq} = $newsp;
          delete $lineage{$clsp};
          $lineage{$newsp} = $lineage{$refsp};
        }
      }
      foreach my $rseqID ( @refseqID)
      {
        {
          $newtx{$rseqID} = $newsp;
        }
      }
    }
      ## If the reference sp has more seq's, then we need to 
      ## 1. update the taxonomy (seqid -> tx) to be refsp_clsp
      ## 2. delete the clsp from the lineage hash
      ## 3. add a new key/value with the refsp_clsp to the ref_lineage
    elsif ($winner == 1)
    {
      $newsp = "$clsp"."_"."$refsp";
      chomp $newsp;
      if ($debug)
      {
        print "$refsp is now $newsp\n";
      }
      foreach my $rseqID ( @refseqID)
      {
        {
          $newtx{$rseqID} = $newsp;
          delete $lineage{$refsp};
          $lineage{$newsp} = $lineage{$clsp};
        }
      }
      foreach my $clseq ( @clSeqs)
      {
        {
          $newtx{$clseq} = $newsp;
        }
      }
    }
  }
  else 
  {
    next;
  }
  %doneSeqs = map { $_ => 1 } @clSeqs;
}

my $fixedTxFile = "$outDir/spp_merged.tx";
print "---Writing merged taxonomy to $fixedTxFile\n";
open OUT, ">$fixedTxFile" or die "Cannot open $fixedTxFile for appending: $OS_ERROR";
for my $keys (keys %fixedtx)
{
  print OUT "$keys\t$fixedtx{$keys}\n";
}
close OUT;

my $newLineageFile = "$outDir/spp_merged.lineage";
print "---Writing merged lineages to $newLineageFile\n";
open OUT, ">$newLineageFile" or die "Cannot open $newLineageFile for appending: $OS_ERROR";
for my $keys (keys %newlineage)
{
  print OUT "$keys\t$newlineage{$keys}\n";
}
close OUT;

exit;