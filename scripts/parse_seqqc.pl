#!/usr/bin/perl
#uploadqc.pl

open OUT, ">sequence.stats.txt" or die $!;
print OUT join("\t",'Sample','total.raw','total.trimmed','ontarget','maprate','propair','Percent.ontarget','Mean.MapQ',
	       'Percent.QualReads','Median.Mismatch','Indel.Rate','Error.Rate','Percent.Dups','Library.Size',
	       'Median.Insert','Mean.Insert','Std.insert',
	       'all.avg.depth','all.median.depth','all.perc.100x','all.perc.200x','all.perc.500x',
	       'dedup.avg.depth','dedup.median.depth','dedup.perc.100x','dedup.perc.200x','dedup.perc.500x'),"\n";

my @statfiles = @ARGV;

foreach $sfile (@statfiles) {
  $sfile =~ m/(\S+)\.genomecov.txt/;
  my $prefix = $1;
  my %hash;
  open FLAG, "<$prefix\.trimreport.txt" or die $!;
  while (my $line = <FLAG>) {
    chomp($line);
    my ($file,$raw,$trim) = split(/\t/,$line);
    $hash{rawct} += $raw;
    $hash{trimct} += $trim;
  }
  open FLAG, "<$prefix\.flagstat.txt" or die $!;
  while (my $line = <FLAG>) {
    chomp($line);
    if ($line =~ m/(\d+) \+ \d+ in total/) {
      $hash{total} = $1;
      $hash{total} = $hash{trimct} if ($hash{trimct} && $hash{trimct} > $hash{total});
    }elsif ($line =~ m/(\d+) \+ \d+ read1/) {
      $hash{pairs} = $1;
    }elsif ($line =~ m/(\d+) \+ \d+ mapped/) {
      $hash{maprate} = 100*sprintf("%.4f",$1/$hash{total});
    }elsif ($line =~ m/(\d+) \+ \d+ properly paired/) {
      $hash{propair} = 100*sprintf("%.4f",$1/$hash{total});
    }elsif ($line =~ m/(\d+) \+ \d+ properly paired/) {
      $hash{propair} = 100*sprintf("%.4f",$1/$hash{total});
	}
  }
  # open FLAG, "<$prefix\.ontarget.flagstat.txt" or die $!;
  # while (my $line = <FLAG>) {
  #   chomp($line);
  #   if ($line =~ m/(\d+) \+ \d+ in total/) {
  #     $hash{ontarget} = $1;
  #   }elsif ($line =~ m/(\d+) \+ \d+ read1/) {
  #     $hash{ontarget_pairs} = $1;
  #   }elsif ($line =~ m/(\d+) \+ \d+ properly paired/) {
  #     $hash{ontarget_propair} = 100*sprintf("%.4f",$1/$hash{total});
  #   }elsif ($line =~ m/(\d+) \+ \d+ properly paired/) {
  #     $hash{ontarget_propair} = 100*sprintf("%.4f",$1/$hash{total});
  #   }
  # }
  open ASM, "<$prefix\.alignmentsummarymetrics.txt" or die $!;
  while (my $line = <ASM>) {
    chomp($line);
    if ($line =~ m/## METRICS/) {
      $header = <ASM>;
      chomp($header);
      my @stats = split(/\t/,$header);
      while (my $nums = <ASM>) {
	chomp($nums);
	next unless ($nums =~ m/^PAIR/);
	my @row = split(/\t/,$nums);
	my %info;
	foreach my $i (0..$#stats) {
	  $info{$stats[$i]} = $row[$i];
	}
	$hash{ontarget} = $info{PF_READS_ALIGNED};
	$hash{qualreads} = $info{PF_HQ_ALIGNED_READS};
	$hash{median_mismatch} = $info{PF_MISMATCH_RATE};
	$hash{indelrate} = $info{PF_INDEL_RATE};
	$hash{error_rate} = $info{PF_HQ_ERROR_RATE};
      }
    }
  }
  close ASM;
  open MM, "<$prefix\.meanmap.txt" or die $!;
  my $meanmap = <MM>;
  chomp($meanmap);
  $meanmap =~ s/Mean MAPQ = //;
  close MM;
  open DUP, "<$prefix\.libcomplex.txt" or die $!;
  while (my $line = <DUP>) {
    chomp($line);
    if ($line =~ m/## METRICS/) {
      $header = <DUP>;
      $nums = <DUP>;
      chomp($header);
      chomp($nums);
      my @stats = split(/\t/,$header);
      my @row = split(/\t/,$nums);
      my %info;
      foreach my $i (0..$#stats) {
	$info{$stats[$i]} = $row[$i];
      }
      $hash{percdups} = sprintf("%.4f",$info{PERCENT_DUPLICATION});
      $hash{libsize} = $info{ESTIMATED_LIBRARY_SIZE};
    }
  }
  close DUP;
  $hash{medinsert} = 0;
  $hash{avginsert} = 0;
  $hash{stdinsert} = 0;
  open DUP, "<$prefix\.hist.txt";
  while (my $line = <DUP>) {
    chomp($line);
    if ($line =~ m/## METRICS/) {
      $header = <DUP>;
      $nums = <DUP>;
      chomp($header);
      chomp($nums);
      my @stats = split(/\t/,$header);
      my @row = split(/\t/,$nums);
      my %info;
      foreach my $i (0..$#stats) {
	$info{$stats[$i]} = $row[$i];
      }
      $hash{medinsert} = sprintf("%.0f",$info{MEDIAN_INSERT_SIZE});
      $hash{avginsert} = sprintf("%.0f",$info{MEAN_INSERT_SIZE});
      $hash{stdinsert} = sprintf("%.0f",$info{STANDARD_DEVIATION});
    }
  }
  my %cov;
  open COV, "<$prefix\.genomecov.txt" or die $!;
  my $sumdepth;
  my $totalbases;
  while (my $line = <COV>) {
    chomp($line);
    my ($all,$depth,$bp,$total,$percent) = split(/\t/,$line);
    $cov{$depth} = $percent;
    $sumdepth += $depth*$bp;
    $totalbases = $total if ($depth == 0);
  }
  my $avgdepth = sprintf("%.0f",$sumdepth/$totalbases);
  my @depths = sort {$a <=> $b} keys %cov;
  my @perc = @cov{@depths};
  my @cum_sum = cumsum(@perc);
  my $median = 0;
  foreach my $i (0..$#cum_sum) {
    if ($cum_sum[$i] < 0.5) {
      $median = $i;
	}
  }
  #$hash{'perc10x'} = 100*sprintf("%.4f",1-$cum_sum[10]);
  #$hash{'perc50x'} = 100*sprintf("%.4f",1-$cum_sum[50]);
  $hash{'perc100x'} = 100*sprintf("%.4f",1-$cum_sum[100]);
  $hash{'perc200x'} = 100*sprintf("%.4f",1-$cum_sum[200]);
  $hash{'perc500x'} = 100*sprintf("%.4f",1-$cum_sum[500]);
  open COV, "<$prefix\.dedupcov.txt" or die $!;
  $sumdepth = 0;
  $totalbases = 0;
  my %dedup_cov;
  while (my $line = <COV>) {
    chomp($line);
    my ($all,$depth,$bp,$total,$percent) = split(/\t/,$line);
    $dedup_cov{$depth} = $percent;
    $sumdepth += $depth*$bp;
    $totalbases = $total if ($depth == 0);
  }
  my $dedup_avgdepth = sprintf("%.0f",$sumdepth/$totalbases);
  my @depths = sort {$a <=> $b} keys %dedup_cov;
  my @perc = @dedup_cov{@depths};
  my @cum_sum = cumsum(@perc);
  my $dedup_median = 0;
  foreach my $i (0..$#cum_sum) {
    if ($cum_sum[$i] < 0.5) {
      $dedup_median = $i;
    }
  }
  #$hash{'dedup.perc10x'} = 100*sprintf("%.4f",1-$cum_sum[10]);
  #$hash{'dedup.perc50x'} = 100*sprintf("%.4f",1-$cum_sum[50]);
  $hash{'dedup.perc100x'} = 100*sprintf("%.4f",1-$cum_sum[100]);
  $hash{'dedup.perc200x'} = 100*sprintf("%.4f",1-$cum_sum[200]);
  $hash{'dedup.perc500x'} = 100*sprintf("%.4f",1-$cum_sum[500]);

  print OUT join("\t",$prefix,$hash{rawct},$hash{total},$hash{ontarget},$hash{maprate},$hash{propair},
		 sprintf("%.2f",100*$hash{ontarget}/$hash{total}),$meanmap,
		 sprintf("%.2f",100*$hash{qualreads}/$hash{ontarget}),$hash{median_mismatch},
		 $hash{indelrate},$hash{error_rate},sprintf("%.2f",100*$hash{percdups}),$hash{libsize},
		 $hash{medinsert},$hash{avginsert},$hash{stdinsert},$avgdepth,$median,$hash{'perc100x'},
		 $hash{'perc200x'},$hash{'perc500x'},$dedup_avgdepth,$dedup_median,
		 $hash{'dedup.perc100x'},$hash{'dedup.perc200x'},$hash{'dedup.perc500x'}),"\n";
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
