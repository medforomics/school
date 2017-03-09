#!/usr/bin/perl -w
#get_regions_coverage.pl

open ANNOT, "</project/shared/bicf_workflow_ref/GRCh38/utswv2_cds.oncomine_cosmic.txt" or die $!;
my %annot;
while (my $line = <ANNOT>) {
  chomp($line);
  my ($chrom,$start,$end,$oncomine,$cosmic_census,
      $cosmic_hotspot) = split(/\t/,$line);
  my $region = $chrom.":".$start."-".$end;
  $annot{$region} = [$oncomine,$cosmic_census,$cosmic_hotspot];
}
close ANNOT;

my $cfile = shift @ARGV;
my @fullpath = split(/\//,$cfile);
my $fname = pop @fullpath;
my $dir = join("/",@fullpath);
$prefix = (split(/\./,$fname))[0];

open OUT2, ">$dir\/$prefix\.cdscovstats.txt" or die $!;
print OUT2 join("\t",'Sample','Chromosome','Start','End','ExonName',
		'Oncomine','CosmicConsesus','CosmicHotspot',
		'MinDepth','MaxDepth','MedianDepth','AvgDepth',
		'Fraction100X+','Fraction20X-','BP100X+','BP20X-','TotalBP'),"\n";

open MISS, ">$dir\/$prefix\.missingcov.txt" or die $!;
print MISS join("\t",'Sample','Chromosome','Start','End','ExonName',
		'Oncomine','CosmicConsesus','CosmicHotspot',
		'MinDepth','MaxDepth','MedianDepth','AvgDepth',
		'Fraction100X+','Fraction20X-','BP100X+','BP20X-','TotalBP'),"\n";

open COV, "<$cfile" or die $!;
while (my $line = <COV>) {
  chomp($line);
  next if ($line =~ m/^all/);
  my ($chrom,$pos,$end,$name,$depth,$bp,$total,$percent) = split(/\t/,$line);
  my $region = $chrom.":".$pos."-".$end;
  $cov{$region}{$depth} = $percent;
  $covstats{$region}{sumdepth} += $depth*$bp;
  $exons{$region} = $name;
  $totalbp{$region} = $total;
  if ($depth <= 20) {
    $poor{$region} += $bp;
  }
  if ($depth >= 100) {
    $good{$region} += $bp;
  }
}
close COV;
foreach $reg (keys %totalbp) {
  $good{$reg} = 0 unless $good{$reg};
  $poor{$reg} = 0 unless $poor{$reg};
  if ($totalbp{$reg} < 1) {
    next;
  }
  my $avgdepth = sprintf("%.0f",$covstats{$reg}{sumdepth}/$totalbp{$reg});
    my @depths = sort {$a <=> $b} keys %{$cov{$reg}};
  %covhash = %{$cov{$reg}};
  my @perc = @covhash{@depths};
  my @cum_sum = cumsum(@perc);
  my $mediandep = 0;
  foreach my $i (0..$#cum_sum) {
    if ($cum_sum[$i] < 0.5) {
      $mediandep = $depths[$i];
    }
  }
  $good_fract = sprintf("%.4f",$good{$reg}/$totalbp{$reg});
  $poor_fract = sprintf("%.4f",$poor{$reg}/$totalbp{$reg});
  $annot{$reg} = [0,0,0] unless $annot{$reg};
  print OUT2 join("\t",$prefix,split(/:|-/,$reg),$exons{$reg},@{$annot{$reg}},$depths[0],$depths[-1],$mediandep,$avgdepth,$good_fract,$poor_fract,$good{$reg},$poor{$reg},$totalbp{$reg}),"\n";
  if ($depths[0] < 100) {
    print MISS join("\t",$prefix,split(/:|-/,$reg),$exons{$reg},@{$annot{$reg}},$depths[0],$depths[-1],$mediandep,$avgdepth,$good_fract,$poor_fract,$good{$reg},$poor{$reg},$totalbp{$reg}),"\n";
  }
}
close MISS;
close OUT2;

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
