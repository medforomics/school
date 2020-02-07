#!/usr/bin/perl -w
#sequenceqc_rnaseq.pl

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %opt = ();
my $results = GetOptions (\%opt,'refdir|r=s','help|h');

my @files = @ARGV;
chomp(@files);

#### Begin Version Information ######
my $gittag = `cd /project/PHG/PHG_Clinical/clinseq_workflows;git describe --abbrev=0 --tags`;
chomp $gittag;
#### End Version Information ######

foreach my $file (@files) {
  chomp($file);
  $file =~ m/(\S+)\.flagstat.txt/;
  my @directory = split(/\//, $file);
  $prefix = $directory[-1];
  my $sample = (split(/\./,$prefix))[0];
  open OUT, ">".$sample.".sequence.stats.txt" or die;#separateFilesPerSample
  $hisatfile = $file;
  $hisatfile =~ s/flagstat/alignerout/;
  open HISAT, "<$hisatfile";
  my ($total, $pairs,$read2ct,$maprate,$concorrate);
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
  
  my @stats=stat($file);
  my ($day,$month,$year) = (localtime($stats[9]))[3,4,5];
  $year += 1900;
  $month++;
  $date = join("-",$year,sprintf("%02s",$month),sprintf("%02s",$day));
  $fileowner = 's'.$stats[4];

  $total = $pairs*2 if ($total == $pairs);
  my $status = 'PASS';
  $status = 'FAIL' if ($maprate < 0.90);
  $status = 'FAIL' if ($total < 6000000);  
  
  print OUT join("\n","Sample\t".$sample,"Total_Raw_Count\t".$total, "Read1_Map\t".$pairs,"Read2_Map\t".$read2ct,
    "Map_Rate\t".$maprate,"Concordant_Rate\t".$concorrate,"Alignment_Status\t".$status,"Alignment_Date\t".$date,
    "File_Owner\t".$fileowner,"Workflow_Version\t".$gittag),"\n";
  system(qq{cat $opt{refdir}\/reference_info.txt >> $prefix\.sequence.stats.txt});
}
