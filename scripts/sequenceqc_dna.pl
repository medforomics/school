#!/usr/bin/perl -w
#sequenceqc_alignment.p

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %opt = ();
my $results = GetOptions (\%opt,'refdir|r=s','help|h','gitdir|e=s','user|u=s');

my @statfiles = @ARGV;

my $fileowner = $opt{user};
my $gittag = 'v5';
if ($opt{gitdir}) {
    $gittag = $opt{gitdir};
}elsif($ENV{'gitv'}) {
    $gittag = $ENV{'gitv'};
}
foreach $sfile (@statfiles) {
  $sfile =~ m/(\S+)\.genomecov.txt/;
  my $prefix = $1;
  my %hash;
  open FLAG, "<$prefix\.trimreport.txt";
  while (my $line = <FLAG>) {
    chomp($line);
    my ($file,$raw,$trim) = split(/\t/,$line);
    $hash{rawct} += $raw;
  }
  open FLAG, "<$prefix\.flagstat.txt" or die $!;
  while (my $line = <FLAG>) {
    chomp($line);
    if ($line =~ m/(\d+) \+ \d+ in total/) {
      $hash{total} = $1;
    }elsif ($line =~ m/(\d+) \+ \d+ read1/) {
      $hash{pairs} = $1;
    }elsif ($line =~ m/(\d+) \+ \d+ mapped/) {
      $hash{maprate} = 100*sprintf("%.4f",$1/$hash{total});
    }elsif ($line =~ m/(\d+) \+ \d+ properly paired/) {
      $hash{propair} = 100*sprintf("%.4f",$1/$hash{total});
    }
  }
  close FLAG;
  open FLAG, "<$prefix\.ontarget.flagstat.txt" or die $!;
  while (my $line = <FLAG>) {
    chomp($line);
    if ($line =~ m/(\d+) \+ \d+ in total/) {
      $hash{ontarget} = $1;
    }
  }
  unless ($hash{rawct}) {
    $hash{rawct} = $hash{total};
  }
  my %lc;
  open DUP, "<$prefix.libcomplex.txt" or die $!;
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
      $lc{TOTREADSLC} += $info{UNPAIRED_READS_EXAMINED} + $info{READ_PAIRS_EXAMINED};
      $hash{libsize} = $info{ESTIMATED_LIBRARY_SIZE};
      $lc{TOTDUPLC} +=  $info{UNPAIRED_READ_DUPLICATES} + $info{READ_PAIR_DUPLICATES};
    }
  }
  close DUP;
  $hash{percdups} = sprintf("%.4f",$lc{TOTDUPLC}/$lc{TOTREADSLC});
  my %cov;
  open COV, "<$prefix\.genomecov.txt" or die $!;
  my $sumdepth;
  my $totalbases;
  while (my $line = <COV>) {
    chomp($line);
    my ($all,$depth,$bp,$total,$percent) = split(/\t/,$line);
    $cov{$depth} = $percent;
    $sumdepth += $depth*$bp;
    $totalbases = $total unless ($totalbases);
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
  
  #### Begin File Information ########
  my $status = 'PASS';
  $status = 'FAIL' if ($hash{maprate} < 0.90 && $hash{'perc100x'} < 0.90);
  my @stats = stat("$prefix\.flagstat.txt");
  my ($day,$month,$year) = (localtime($stats[9]))[3,4,5];
  $year += 1900;
  $month++;
  $date = join("-",$year,sprintf("%02s",$month),sprintf("%02s",$day));
  $fileowner = 's'.$stats[4] unless $fileowner;
  $hash{status} = $status;
  $hash{date}=$date;
  $hash{fileowner} = $fileowner;
  ##### End File Information ########
  ##### START separateFilesPerSample ######
  open OUT, ">".$prefix.".sequence.stats.txt" or die $!;
  print OUT join("\n", "Sample\t".$prefix,"Total_Raw_Count\t".$hash{rawct},
		 "On_Target\t".$hash{ontarget},"Map_Rate\t".$hash{maprate},
		 "Properly_Paired\t".$hash{propair},
		 "Percent_Duplicates\t".sprintf("%.2f",100*$hash{percdups}),
		 "Percent_on_Target\t".sprintf("%.2f",100*$hash{ontarget}/$hash{rawct}),
		 "All_Average_Depth\t".$avgdepth,"All_Median_Depth\t".$median,
		 "Percent_over_100x\t".$hash{'perc100x'},
		 "Percent_over_200x\t".$hash{'perc200x'},
		 "Percent_over_500x\t".$hash{'perc500x'},
		 "Alignment_Status\t".$hash{status},"Alignment_Date\t".$hash{date},
		 "File_Owner\t".$hash{fileowner},"Workflow_Version\t".$gittag),"\n";
  close OUT;
  system(qq{cat $opt{refdir}\/reference_info.txt >> $prefix\.sequence.stats.txt});
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
