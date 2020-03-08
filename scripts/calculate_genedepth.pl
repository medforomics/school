#!/usr/bin/perl -w
#get_regions_coverage.pl

my $cfile = shift @ARGV;
my $prefix = (split(/\./,$cfile))[0];

open OUT, ">$prefix\.genecoverage.txt" or die $!;
print OUT join("\t",'Sample','Gene','MinDepth','MaxDepth',
	       'MedianDepth','AvgDepth','TotalBP'),"\n";
my %cov;
my %covstats;
my %exons;
my %totalbp;
my %good;

open COV, "<$cfile" or die $!;
while (my $line = <COV>) {
  chomp($line);
  next if ($line =~ m/^all/);
  my ($chrom,$pos,$end,$name,$depth,$bp,$total,$percent) = split(/\t/,$line);
  my ($gene,$exonname) = split(/:/,$name);
  $cov{$gene}{$depth} += $bp;
  $covstats{$gene}{sumdepth} += $depth*$bp;
  $totalbp{$gene}{$name} = $total;
  if ($depth >= 100) {
    $good{$gene} += $bp;
  }
}
close COV;
foreach $reg (keys %totalbp) {
  $good{$reg} = 0 unless $good{$reg};
  my $totbp = 0;
  foreach $exon (keys %{$totalbp{$reg}}) {
    $totbp += $totalbp{$reg}{$exon};
  }
  if ($totbp < 1) {
    next;
  }
  my $avgdepth = sprintf("%.0f",$covstats{$reg}{sumdepth}/$totbp);
  my @depths = sort {$a <=> $b} keys %{$cov{$reg}};
  #my %covhash = %{$cov{$reg}};
  foreach $d (@depths) {
    $covhash{$d} = $cov{$reg}{$d}/$totbp;
  }
  my @perc = @covhash{@depths};
  my @cum_sum = cumsum(@perc);
  my $mediandep = 0;
  foreach my $i (0..$#cum_sum) {
    if ($cum_sum[$i] < 0.5) {
      $mediandep = $depths[$i];
    }
  }
  #$good_fract = sprintf("%.4f",$good{$reg}/$totbp);
  print OUT join("\t",$prefix,split(/:/,$reg),$depths[0],$depths[-1],
		 $mediandep,$avgdepth,$totbp),"\n";
}


sub cumsum {
  my @nums = @_;
  my @cumsum = ();
  my $mid = 0;
  for my $i (0..$#nums) {
    $mid += $nums[$i];
    push(@cumsum,$mid);
  }
  return @cumsum;
}
