#!/usr/bin/perl -w
#parse_alignment_stats.pl

my @files = @ARGV;
chomp(@files);

my $gittag = `cd /project/PHG/PHG_Clinical/clinseq_workflows;git describe --tag`;
chomp $gittag;

foreach my $file (@files) {
  chomp($file);
  $file =~ m/(\S+)\.flagstat.txt/;
  my @directory = split(/\//, $file);
  $prefix = $directory[-1];
  my $sample = (split(/\./,$prefix))[0];
  open OUT, ">".$sample."_sequence.stats.txt" or die;#separateFilesPerSample
  $hisatfile = $file;
  $hisatfile =~ s/flagstat/hisatout/;
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
  #### Begin File Information ########
  my $status = 'PASS';
  $status = 'FAIL' if ($maprate < 0.90);
  my @stats=stat($file);
  my ($day,$month,$year) = (localtime($stats[9]))[3,4,5];
  $year += 1900;
  $month++;
  $date = join("-",$year,sprintf("%02s",$month),sprintf("%02s",$day));
  $fileowner = 's'.$stats[4];
 ##### End File Information ########
 ##### Begin Reference Information ########

  $filtvcffile = $file;
  $filtvcffile =~ s/flagstat\.txt/filt\.vcf/;
  open FILT, "<$filtvcffile";
   my ($cosmic_ref,$dbsnp_ref,$gen_ref)= ("") x 3;
   while (my $line = <FILT>) {
    chomp($line);
    if ($line =~ m/cosmic*/) {
  $line =~ s/^##.*\s([\S]+cosmic[\S]+).*/$1/;
        $cosmic_ref = $line;
    }
    if ($line =~ m/source=.*dbSnp/) {
  $line =~ s/^.*source=([\S]+dbSnp[\S]+).*/$1/;
  $line =~ s/\)//;
        $dbsnp_ref = $line;
    }
    if ($line =~ m/##reference=file:\/\//) {
  $line =~ s/##reference=file:\/\/([\S]+)*/$1/;
        $gen_ref = $line;
    }

  }
 ##### End Reference Information ########

  $total = $pairs*2 if ($total == $pairs);
  print OUT join("\n","Sample\t".$sample,"Total_Raw_Count\t".$total, "Read1_Map\t".$pairs,"Read2_Map\t".$read2ct,
    "Map_Rate\t".$maprate,"Concordant_Rate\t".$concorrate,"Alignment_Status\t".$status,"Alignment_Date\t".$date,
    "File_Owner\t".$fileowner,"Workflow_Version\t".$gittag,"Cosmic_Reference\t".$cosmic_ref,
    "dbSnp_Reference\t".$dbsnp_ref, "Genome_Reference\t".$gen_ref),"\n";
}
