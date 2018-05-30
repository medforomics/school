#!/usr/bin/perl -w
#integrate_datasets.pl

#module load vcftools/0.1.14 samtools/1.6 bedtools/2.26.0 
my ($subject,$tumorid,$normalid,$refdir,$rnaseq_vcf,$rnaseq_ntct) = @ARGV;

open OM, "<$refdir\/clinseq_prj/panel1385.genelist.txt" or die $!;
while (my $line = <OM>) {
  chomp($line);
  $keep{$line} = 1;
}
close OM;
open OM, "<$refdir\/clinseq_prj/cancer.genelist.txt" or die $!;
while (my $line = <OM>) {
  chomp($line);
  $cgenelist{$line} = 1;
}
close OM;

my $rnaseqid;
if ($rnaseq_ntct && -e $rnaseq_ntct) {
  open NRC, "<$rnaseq_ntct" or die $!;
  while (my $line = <NRC>) {
    chomp($line);
    my ($chrom,$pos,$ref,$depth,@reads) = split(/\t/,$line);
    next unless ($depth > 10);
    $chrom = 'chr'.$chrom if ($chrom !~ m/^chr/);
    my $ro;
    my %hash;
    foreach my $rct (@reads) {
	my ($base,$basect,@otherstats) = split(/:/,$rct);
	if ($ref eq $base) {
	    $hash{$base} = $basect;
	}else {
	  if ($base =~ m/\+|\-/) {
	    $base =~ s/\+/$ref/;
	    #$base =~ s/\-/$ref/;
	  }
	  $hash{$base} = $basect if ($basect);
	}
      }
    $rnaval{$chrom}{$pos} = [\%hash,$depth];
  }
  close NRC;
  open RVCF, "gunzip -c $rnaseq_vcf |" or die $!;
 W1:while (my $line = <RVCF>) {
    chomp($line);
    if ($line =~ m/^#CHROM/) {
      my @header = split(/\t/,$line);
      ($chrom, $pos,$id,$ref,$alt,$score,
       $filter,$info,$format,@gtheader) = split(/\t/, $line);
      $rnaseqid = $gtheader[0];
    }
    if ($line =~ m/^#/) {
      next;
    }
    my ($chrom, $pos,$id,$ref,$alt,$score,
	$filter,$annot,$format,@gts) = split(/\t/, $line);
    next if ($ref =~ m/\./ || $alt =~ m/\./ || $alt=~ m/,X/);
    $chrom = 'chr'.$chrom if ($chrom !~ m/^chr/);
    my %hash = ();
    foreach $a (split(/;/,$annot)) {
	my ($key,$val) = split(/=/,$a);
	$hash{$key} = $val unless ($hash{$key});
    }
    my @deschead = split(/:/,$format);
    my $allele_info = shift @gts;
    @ainfo = split(/:/, $allele_info);
    my %gtinfo = ();
    my @mutallfreq = ();
    foreach my $k (0..$#ainfo) {
	$gtinfo{$deschead[$k]} = $ainfo[$k];
    }
    $gtinfo{DP} = (split(/,/,$gtinfo{DP}))[0];
    next W1 if ($gtinfo{DP} < 10);
    my ($ro,@altct) = split(/,/,$gtinfo{AD});
    my @alts = split(/,/,$alt);
    my %allct;
    foreach my $j (0..$#altct) {
	$act = $altct[$j];
	$base = $alts[$j];
      $allct{$base} = $act;
    }
    $rnaval{$chrom}{$pos} = [\%allct,$gtinfo{DP}];
  }
}

open OUT, ">$subject\.all.vcf" or die $!;
open PASS, ">$subject\.pass.vcf" or die $!;

my @sampids;
open IN, "gunzip -c somatic_germline.vcf.gz |" or die $!;
W1:while (my $line = <IN>) {
  chomp($line);
  if ($line =~ m/^#CHROM/) {
    print OUT qq{##INFO=<ID=RnaSeqAF,Number=A,Type=Float,Description="RNASeq Allele Frequency">\n};
    print OUT qq{##INFO=<ID=RnaSeqDP,Number=1,Type=Integer,Description="RNASeq read depth">\n};
    print PASS qq{##INFO=<ID=RnaSeqAF,Number=A,Type=Float,Description="RNASeq Allele Frequency">\n};
    print PASS qq{##INFO=<ID=RnaSeqDP,Number=1,Type=Integer,Description="RNASeq read depth">\n};
    my @header = split(/\t/,$line);
    ($chrom, $pos,$id,$ref,$alt,$score,
     $filter,$info,$format,@gtheader) = split(/\t/, $line);
    @sampids = ($tumorid,$normalid);
    push @sampids, $rnaseqid if ($rnaseqid);
    print OUT join("\t",$chrom,$pos,$id,$ref,$alt,$score,
		   $filter,$info,$format,@sampids),"\n";
    print PASS join("\t",$chrom,$pos,$id,$ref,$alt,$score,
		    $filter,$info,$format,@sampids),"\n";
    next;
  } elsif ($line =~ m/^#/) {
    print OUT $line,"\n";
    print PASS $line,"\n";
    next;
  }
  my ($chrom, $pos,$id,$ref,$alt,$score,
      $filter,$annot,$format,@gts) = split(/\t/, $line);
  next if ($ref =~ m/\./ || $alt =~ m/\./ || $alt=~ m/,X/);
  my %hash = ();
  foreach $a (split(/;/,$annot)) {
    my ($key,$val) = split(/=/,$a);
    $hash{$key} = $val unless ($hash{$key});
  }
  $hash{'HG38Loci'} = join(":",$chrom,$pos);
  my %fail;
  $fail{'UTSWBlacklist'} = 1 if ($hash{UTSWBlacklist});
  my @exacaf;
  my $exacaf;
  if ($hash{AF_POPMAX}) {
    foreach (split(/,/,$hash{AF_POPMAX})) {
      push @exacaf, $_ if ($_ ne '.');
    }
    @exacaf = sort {$b <=> $a} @exacaf;
    $exacaf = $exacaf[0] if ($exacaf[0]);
  }elsif ($hash{AC_POPMAX} && $hash{AN_POPMAX}) {
    my @exacs = split(/,/,$hash{AC_POPMAX});
    my $ac = 0;
    foreach $val (@exacs) {
      $ac += $val if ($val =~ m/^\d+$/);
    }
    my @exans = split(/,/,$hash{AN_POPMAX});
    my $an = 0;
    foreach $val (@exans) {
      $an += $val if ($val =~ m/^\d+$/);
    }
    $exacaf = sprintf("%.4f",$ac/$an) if ($ac > 0 && $an > 10);
  }
  $fail{'COMMON'} = 1 if ($exacaf && $exacaf > 0.01);
  $fail{'StrandBias'} = 1 if (($hash{FS} && $hash{FS} > 60) || $filter =~ m/strandBias/i);
  my $cosmicsubj = 0;
  if ($hash{CNT}) {
    my @cosmicct = split(/,/,$hash{CNT}); 
    foreach $val (@cosmicct) {
      $cosmicsubj += $val if ($val =~ m/^\d+$/);
    }
  }
  my %gtinfo = ();
  my @deschead = split(/:/,$format);
 F1:foreach my $k (0..$#gtheader) {
    my $subjid = $gtheader[$k];
    my $allele_info = $gts[$k];
    my @ainfo = split(/:/, $allele_info);
    my @mutallfreq = ();
    foreach my $k (0..$#ainfo) {
      $gtinfo{$subjid}{$deschead[$k]} = $ainfo[$k];
      $hash{$deschead[$k]} = $ainfo[$k] if ($subjid eq $tumorid);
    }
    $gtinfo{$subjid}{DP} = (split(/,/,$gtinfo{$subjid}{DP}))[0] if ($gtinfo{$subjid}{DP});
    next F1 unless ($gtinfo{$subjid}{DP} && $gtinfo{$subjid}{DP} ne '.' && $gtinfo{$subjid}{DP} >= 1);
    my @altct = split(/,/,$gtinfo{$subjid}{AD});
    my $refct = shift @altct;
    @altct2 = split(/,/,$gtinfo{$subjid}{AO});
    if (scalar(@altct) ne scalar(@altct2)) {
	warn "Inconsistent Allele counts @ $chrom,$pos,$alt,$gtinfo{$subjid}{AD},$gtinfo{$subjid}{AO}\n";
    }
    my $total = $refct;
    foreach  my $act (@altct) {
	next if ($act eq '.');
	$total += $act;
	push @mutallfreq, sprintf("%.4f",$act/$gtinfo{$subjid}{DP});
    }
    unless ($gtinfo{$subjid}{DP} eq $total) {
	warn "Inconsistent Depth counts @ $chrom,$pos,$alt,$gtinfo{$subjid}{AD},$gtinfo{$subjid}{DP}\n";
    }
    $gtinfo{$subjid}{MAF} = \@mutallfreq;
  }
  next unless ($gtinfo{$tumorid}{DP} && $gtinfo{$tumorid}{DP} ne '.' && $gtinfo{$tumorid}{DP} >= 20);
  unless ($gtinfo{$tumorid}{AO} =~ m/\d+/ && $gtinfo{$tumorid}{AD} =~ m/,/) {
      warn "Missing Alt:$line\n";
  }
  @tumormaf = @{$gtinfo{$tumorid}{MAF}};
  @tumoraltct = split(/,/,$gtinfo{$tumorid}{AO});
  
  if (exists $hash{INDEL}) {
      $hash{TYPE} = 'indel';
  }
  $hash{TYPE} = 'ambi' unless ($hash{"TYPE"});
  next if ($tumoraltct[0] eq '.');
  $hash{AF} = join(",",@tumormaf);
  my @callers ;
  @callers = split(/,/, $hash{CallSet}) if ($hash{CallSet});
  if ($hash{SomaticCallSet}) {
    @callers = (@callers,split(/,/,$hash{SomaticCallSet}));
  }
  $hash{CallSet} = join(",",@callers);
  if ($id =~ m/COS/ && $cosmicsubj >= 5) {
    $fail{'LowAltCt'} = 1 if ($tumoraltct[0] < 3);
    $fail{'LowMAF'} = 1 if ($tumormaf[0] < 0.01);
    #$fail{'LowMAF'} = 1 if ($tumormaf[0] < 0.03 && $hash{TYPE} ne 'snp');
  }else {
    $fail{'OneCaller'} = 1 if (scalar(@callers) < 2);
    $fail{'LowAltCt'} = 1 if ($tumoraltct[0] < 8);
    $fail{'LowMAF'} = 1 if ($tumormaf[0] < 0.05);
    #$fail{'LowMAF'} = 1 if ($tumormaf[0] < 0.10 && $hash{TYPE} ne 'snp');
  }
  delete $hash{SOMATIC};
  $hash{SS} = 5  unless ($hash{SS});
  if ($gtinfo{$normalid} && exists $gtinfo{$normalid}{MAF}) {
      @normalmaf = @{$gtinfo{$normalid}{MAF}};
      $hash{NormalAF} = $normalmaf[0];
      $hash{NormalDP} = $gtinfo{$normalid}{DP};
    if ($normalmaf[0] >= 0.25) {
      $hash{SS} = 1;
    }elsif ($tumormaf[0] < 0.05 && ($normalmaf[0] > 0.01 || $normalmaf[0]*5 > $tumormaf[0])) {
      next;
    }elsif ($normalmaf[0] >= 0.05 || $normalmaf[0]*5 > $tumormaf[0]) {
      $hash{'HighFreqNormalAF'} = 1;
    }else {
      $hash{SS} = 2;
    }
  }
  $gtinfo{$normalid} ={GT=>'.',DP=>'.',AO=>'.',AD=>'.',RO=>'.'} unless ($gtinfo{$normalid});
  if ($rnaval{$chrom}{$pos}) {
    $gtinfo{$rnaseqid} ={GT=>'.',DP=>'.',AO=>'.',AD=>'.',RO=>'.'};
    my ($rnahashref,$rnadp) = @{$rnaval{$chrom}{$pos}};
    if ($rnadp > 10) {
      my %rnantct = %{$rnahashref};
      my @altcts;
      my $totalaltct =0;
      foreach $altnt (split(/,/,$alt)) {
	my $ct = $rnantct{$altnt};
	$ct = 0 unless ($ct);
	$totalaltct += $ct;
	      push @altcts, $ct;
      }
      $hash{RnaSeqDP} = $rnadp;
      $hash{RnaSeqAF} = sprintf("%.4f",$altcts[0]/$rnadp);
      $gtinfo{$rnaseqid}{RO} = $rnadp - $totalaltct;
      $gtinfo{$rnaseqid}{AO} = join(",",@altcts);
      $gtinfo{$rnaseqid}{GT} = '.';
      $gtinfo{$rnaseqid}{DP} = $rnadp;
      $gtinfo{$rnaseqid}{AD} = join(",",$gtinfo{$rnaseqid}{RO},@altcts);
    }
  }
  if ($rnaseqid) {
      if ($gtinfo{$rnaseqid}{AO} && $gtinfo{$rnaseqid}{AO} =~ m/^\d+/ &&
	(split(/,/,$gtinfo{$rnaseqid}{AO}))[0] > 2) {
      $hash{RnaSeqValidation} = 1
    }
  }
  my $newformat = 'GT:DP:AD:AO:RO';
  my @newgt;
  foreach $sample (@sampids) {
    my @gtdata;
    foreach $gt (split(/:/,$newformat)) {
      $gtinfo{$sample}{$gt} = '.' unless (exists $gtinfo{$sample}{$gt});
      push @gtdata, $gtinfo{$sample}{$gt};
    }
    push @newgt, join(":",@gtdata);
  }
  next unless ($hash{ANN});
  my $cancergene = 0;
  my $keepforvcf;
  foreach $trx (split(/,/,$hash{ANN})) {
    my ($allele,$effect,$impact,$gene,$geneid,$feature,
	$featureid,$biotype,$rank,$codon,$aa,$pos_dna,$len_cdna,
	$cds_pos,$cds_len,$aapos,$aalen,$distance,$err) = split(/\|/,$trx);
    next unless ($impact =~ m/HIGH|MODERATE/ || $effect =~ /splice/i);
    next unless $keep{$gene};
    $keepforvcf = $gene;
    $cancergene = 1 if ($cgenelist{$gene});
  }
  next unless $keepforvcf;
  my @fail = sort {$a cmp $b} keys %fail;
  if (scalar(@fail) == 0) {
    $filter = 'PASS';
  }else {
    $filter = join(";", 'FailedQC',@fail);
  }
  my @nannot;
  foreach $info (sort {$a cmp $b} keys %hash) {
    if (defined $hash{$info}) {
      push @nannot, $info."=".$hash{$info};
    }else {
      push @nannot, $info;
    }
  }
  $newannot = join(";",@nannot);
  print PASS join("\t",$chrom, $pos,$id,$ref,$alt,$score,$filter,$newannot,
		  $newformat,@newgt),"\n" if ($filter eq 'PASS');
  print OUT join("\t",$chrom, $pos,$id,$ref,$alt,$score,$filter,$newannot,
		 $newformat,@newgt),"\n" if ($filter eq 'PASS' || $id =~ m/COS/ || $cancergene || $filter eq 'FailedQC;COMMON');
}
close IN;
