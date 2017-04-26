#!/usr/bin/perl -w
#get_regions_coverage.pl

my $cfile = shift @ARGV;
my $prefix = (split(/\./,$cfile))[0];

open OUT2, ">$prefix\_exoncoverage.txt" or die $!;
open MISS, ">$prefix\_lowcoverage.txt" or die $!;
print MISS join("\t",'Sample','Chromosome','Position','ExonName',
		'MinDepth','MaxDepth','MedianDepth','AvgDepth',
		'Fraction100X+','Fraction10X-','BP100X+','TotalBP'),"\n";
print OUT2 join("\t",'Sample','Chromosome','Position','ExonName',
		'MinDepth','MaxDepth','MedianDepth','AvgDepth',
		'Fraction100X+','Fraction10X-','BP100X+','TotalBP'),"\n";
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
  my $region = join(":",$chrom,$pos,$end);
  $cov{$region}{$depth} = $percent;
  $covstats{$region}{sumdepth} += $depth*$bp;
  $exons{$region} = $name;
  $totalbp{$region} = $total;
  if ($depth >= 100) {
    $good{$region} += $bp;
  }
}
close COV;
foreach $reg (keys %totalbp) {
  $good{$reg} = 0 unless $good{$reg};
  if ($totalbp{$reg} < 1) {
    next;
  }
  my $avgdepth = sprintf("%.0f",$covstats{$reg}{sumdepth}/$totalbp{$reg});
  my @depths = sort {$a <=> $b} keys %{$cov{$reg}};
  my %covhash = %{$cov{$reg}};
  my @perc = @covhash{@depths};
  my @cum_sum = cumsum(@perc);
  my $mediandep = 0;
  foreach my $i (0..$#cum_sum) {
    if ($cum_sum[$i] < 0.5) {
      $mediandep = $depths[$i];
    }
  }
  $good_fract = sprintf("%.4f",$good{$reg}/$totalbp{$reg});
  print OUT2 join("\t",$prefix,split(/:/,$reg),$exons{$reg},$depths[0],$depths[-1],$mediandep,$avgdepth,$good_fract,$good{$reg},$totalbp{$reg}),"\n";
  if ($good_fract < 1) {
    print MISS join("\t",$prefix,split(/:/,$reg),$exons{$reg},$depths[0],$depths[-1],$mediandep,$avgdepth,$good_fract,$good{$reg},$totalbp{$reg}),"\n";
  }
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
