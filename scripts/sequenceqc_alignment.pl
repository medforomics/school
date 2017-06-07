#!/usr/bin/perl -w
#uploadqc.pl

my @statfiles = @ARGV;

#### Get version ######
my $gittag = `cd /project/PHG/PHG_Clinical/clinseq_workflows;git describe --tag`;

chomp $gittag;
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
  @depths = sort {$a <=> $b} keys %dedup_cov;
  @perc = @dedup_cov{@depths};
  @cum_sum = cumsum(@perc);
  my $dedup_median = 0;
  foreach my $i (0..$#cum_sum) {
    if ($cum_sum[$i] < 0.5) {
      $dedup_median = $i;
    }
  }
  
  $hash{'dedup.perc100x'} = 100*sprintf("%.4f",1-$cum_sum[100]);
  $hash{'dedup.perc200x'} = 100*sprintf("%.4f",1-$cum_sum[200]);
  $hash{'dedup.perc500x'} = 100*sprintf("%.4f",1-$cum_sum[500]);
  
  #### Begin File Information ########
  my $status = 'PASS';
  $status = 'FAIL' if ($hash{maprate} < 0.90 && $hash{'dedup.perc100x'} < 0.98);
  my @stats = stat("$prefix\.flagstat.txt");
  foreach my $line(@stats){print $line."\n";} 
  my ($day,$month,$year) = (localtime($stats[9]))[3,4,5];
  $year += 1900;
  $month++;
  $date = join("-",$year,sprintf("%02s",$month),sprintf("%02s",$day));
  $fileowner = 's'.$stats[4];
  $hash{status} = $status;
  $hash{date}=$date;
  $hash{fileowner} = $fileowner;
  ##### End File Information ########
  ##### START separateFilesPerSample ######
  open OUT, ">".$prefix."_sequence.stats.txt" or die $!;
  print OUT join("\n", "Sample\t".$prefix,"Total_Raw_Count\t".$hash{rawct},
		 "Total_Trimmed\t".$hash{total},"On_Target\t".$hash{ontarget},"Map_Rate\t".$hash{maprate},
		 "Properly_Paired\t".$hash{propair},"Percent_on_Target\t".sprintf("%.2f",100*$hash{ontarget}/$hash{total}),
		 "Mean_MAPQ\t".$meanmap,"Percent_QualReads\t".sprintf("%.2f",100*$hash{qualreads}/$hash{ontarget}),
		 "Median_Mismatch \t".$hash{median_mismatch},"Indel_Rate\t".$hash{indelrate},
		 "Error_Rate\t".$hash{error_rate},"Percent_Duplicates\t".sprintf("%.2f",100*$hash{percdups}),
		 "Median_Insert_Size\t".$hash{medinsert}, "Average_Insert_Size\t".$hash{avginsert}, 
		 "Std_Insert\t".$hash{stdinsert},"All_Average_Depth\t".$avgdepth, "All_Median_Depth\t".$median,
		 "Percent_over_100x\t".$hash{'perc100x'},"Percent_over_200x\t".$hash{'perc200x'},"Percent_over_500x\t".$hash{'perc500x'},
		 "Dedup_Average_Depth\t".$dedup_avgdepth, "Dedup_Median_Depth\t".$dedup_median,
		 "Dedup_Percent_over_100x\t".$hash{'dedup.perc100x'},"Dedup_Percent_over_200x\t".$hash{'dedup.perc200x'},
		 "Dedup_Percent_over_500x\t".$hash{'dedup.perc500x'},"Alignment_Status\t".$hash{status},"Alignment_Date\t".$hash{date},
		 "File_Owner\t".$hash{fileowner},"Workflow Version\t".$gittag),"\n";
  close OUT;
  ##### END separateFilesPerSample ######
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
