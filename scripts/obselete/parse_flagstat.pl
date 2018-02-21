#!/usr/bin/perl -w
#parse_alignment_stats.pl


my @files = @ARGV;
open OUT, ">alignment.summary.txt" or die $!;
print OUT join("\t",'Sample','Total','R1.map','R2.map','Mapping.Rate','Concordant.Rate'),"\n";
chomp(@files);
foreach my $file (@files) {
  chomp($file);
  my @directory = split(/\//, $file);
  $prefix = $directory[-1];
  my $sample = (split(/\./,$prefix))[0];
  $hisatfile = $file;
  $hisatfile =~ s/flagstat/hisatout/;
  open HISAT, "<$hisatfile";
  my ($total, $pairs,$read1ct,$read2ct,$maprate,$concorrate);
  while (my $line = <HISAT>) {
    chomp($line);
    if ($line =~ m/(\d+) reads; of these:/) {
      $total = $1;
    }elsif ($line =~ m/Number of input reads\s+\|\s+(\d+)/) {
	$total = $1;
    }elsif ($line =~ m/(\d+) \S+ were paired; of these:/) {
      $pairs = $1;
      $total = $pairs*2 if ($total == $pairs);
    }elsif ($line =~ m/(\d+) pairs aligned concordantly 0 times; of these:/) {
      $concorrate = sprintf("%.4f",1-($1/$pairs));
    }elsif ($line =~ m/(\S+)% overall alignment rate/) {
      $maprate = $1/100;
    }
  }
  close HISAT;
  open IN, "<$file" or die $!;
  while ($line = <IN>) {
    chomp($line);
    if ($line =~ m/(\d+) \+ \d+ in total/) {
      $total = $1 unless $total;
    }elsif ($line =~ m/(\d+) \+ \d+ read1/) {
      $pairs = $1;
    }elsif ($line =~ m/(\d+) \+ \d+ read2/) {
      $read2ct = $1;
    }elsif ($line =~ m/(\d+) \+ \d+ mapped\s+\((\S+)%.+\)/) {
      $maprate = sprintf("%.2f",$2/100) unless ($maprate);
    }elsif ($line =~ m/(\d+) \+ \d+ properly paired\s+\((\S+)%.+\)/) {
      $concorrate = sprintf("%.2f",$2/100) unless ($concorrate);
    }
  }
  close IN;
  $total = $pairs*2 if ($total == $pairs);
  print OUT join("\t",$sample,$total, $pairs,$read2ct,$maprate,$concorrate),"\n";
}
