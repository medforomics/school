#!/usr/bin/perl -w
#integrate_datasets.pl

my $refdir = '/project/shared/bicf_workflow_ref/GRCh38/';

open OM, "</project/shared/bicf_workflow_ref/GRCh38/validation.genelist.txt" or die $!;
while (my $line = <OM>) {
  chomp($line);
  $keep{$line} = 1;
}
my ($subject,$samplename,$tumorid,$somatic,$rnaseqid) = @ARGV;
my $inputdir = "/project/PHG/PHG_Clinical/complete/$subject";
system("tabix -f $inputdir\/$tumorid/$tumorid\.annot.vcf.gz");
system("zcat $inputdir\/$tumorid/$tumorid\.annot.vcf.gz > tumor.vcf");
system("/project/PHG/PHG_Clinical/clinseq_workflows/scripts/vcf2bed.pl tumor.vcf |cut -f 1,2,3 > tumor.bed");

if ($rnaseqid ne 'no_rnaseq') {
  system("cat tumor.bed |perl -p -e 's/chr//g' > tumor.nochr.bed");
  system("/project/shared/bicf_workflow_ref/seqprg/bam-readcount/bin/bam-readcount -w 0 -q 10 -b 25 -l tumor.nochr.bed -f /project/shared/bicf_workflow_ref/GRCh38/hisat_genome.fa $inputdir\/$rnaseqid\/$rnaseqid\.final.bam > rnaseq.bamreadct.txt");
  open NRC, "<rnaseq.bamreadct.txt" or die $!;
  while (my $line = <NRC>) {
    chomp($line);
    my ($chr,$pos,$ref,$depth,@reads) = split(/\t/,$line);
    next unless ($depth > 10);
    my $ro;
    my %hash;
    foreach my $rct (@reads) {
      my ($base,$basect,@otherstats) = split(/:/,$rct);
      if ($ref eq $base) {
	$ro = $basect;
      }else {
	$hash{$base} = $basect if ($basect);
      }
    }
    $rnaval{'chr'.$chr}{$pos} = ['','',$depth];
  }
  open IN, "gunzip -c $inputdir\/$rnaseqid\/$rnaseqid\.annot.vcf.gz |" or die $!;
 W1:while (my $line = <IN>) {
    chomp($line);
    if ($line =~ m/^#CHROM/) {
      my @header = split(/\t/,$line);
      ($chrom, $pos,$id,$ref,$alt,$score,
       $filter,$info,$format,@gtheader) = split(/\t/, $line);
      $chrom = 'chr'.$chrom;
    }
    if ($line =~ m/^#/) {
      next;
    }
    my ($chr, $pos,$id,$ref,$alt,$score,
	$filter,$annot,$format,@gts) = split(/\t/, $line);
    my $chrom = 'chr'.$chr;
    next if ($ref =~ m/\./ || $alt =~ m/\./ || $alt=~ m/,X/);
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
    next W1 if ($gtinfo{DP} < 10);
    my ($ro,@altct) = split(/,/,$gtinfo{AD});
    foreach  my $act (@altct) {
      push @mutallfreq, sprintf("%.4f",$act/$gtinfo{DP});
    }
    $maf = $mutallfreq[0];
    $rnaval{$chrom}{$pos} = [$gtinfo{AD},$maf,$gtinfo{DP}];
  }
}
open OUT, ">$inputdir\/$tumorid\.final.vcf" or die $!;
if ($somatic ne 'no_normal') {
  $somatic =~ m/$tumorid\_(.+)/;
  $normal = $1;
  system("/project/shared/bicf_workflow_ref/seqprg/bam-readcount/bin/bam-readcount -w 0 -q 10 -b 25 -l tumor.bed -f /project/shared/bicf_workflow_ref/GRCh38/genome.fa $inputdir\/$normal/$normal\.final.bam > normal.readcts.txt");
  open NRC, "<normal.readcts.txt" or die $!;
  while (my $line = <NRC>) {
    chomp($line);
    my ($chr,$pos,$ref,$depth,@reads)= split(/\t/,$line);
    next if ($depth < 10);
    my $ro;
    my %hash;
    foreach my $rct (@reads) {
      my ($base,$basect,@otherstats) = split(/:/,$rct);
      $hash{$base} = [$basect,$depth,sprintf("%.4f",$basect/$depth)] if ($basect);
    }
    $normals{$chr}{$pos}{NormalDP} = $depth;
    $normals{$chr}{$pos}{BaseCts} = \%hash;
  }
  open NORMVCF, "gunzip -c $inputdir\/$normal/$normal\.annot.vcf.gz |" or die $!;
 W1:while (my $line = <NORMVCF>) {
    chomp($line);
    if ($line =~ m/^#CHROM/) {
      my @header = split(/\t/,$line);
      ($chrom, $pos,$id,$ref,$alt,$score,
       $filter,$info,$format,@gtheader) = split(/\t/, $line);
      next;
    }elsif ($line =~ m/^#/) {
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
    my $allele_info = $gts[0];
    @ainfo = split(/:/, $allele_info);
    my @deschead = split(/:/,$format);
    foreach my $k (0..$#ainfo) {
      $hash{$deschead[$k]} = $ainfo[$k];
    }
    my ($refbasect,@altbasects) = split(/,/,$hash{AD});
    if ($hash{DP} =~ m/,/) {
      $hash{DP} = (split(/,/,$hash{DP}))[0];
    }
    next unless ($hash{DP} && $hash{DP} > 5);
    foreach $altbase (split(/,/,$alt)) {
      my $basect = shift @altbasects;
      next unless ($basect);
      $normals{$chrom}{$pos}{BaseCts}{$altbase} = [$basect,$hash{DP},sprintf("%.4f",$basect/$hash{DP})];
    }
    $normals{$chrom}{$pos}{NormalDP} = $hash{DP};
  }
  open IN, "gunzip -c $inputdir\/$somatic/$somatic\.annot.vcf.gz |" or die $!;
 W1:while (my $line = <IN>) {
    chomp($line);
    if ($line =~ m/^#CHROM/) {
      my @header = split(/\t/,$line);
      ($chrom, $pos,$id,$ref,$alt,$score,
       $filter,$info,$format,@gtheader) = split(/\t/, $line);
      print OUT join("\t",$chrom,$pos,$id,$ref,$alt,$score,$filter,$info,$format,$tumorid),"\n";
      
   } elsif ($line =~ m/^#/) {
      print OUT $line,"\n";
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
    $somline{$chrom}{$pos} = 1;
    my %fail;
    $fail{'UTSWBlacklist'} = 1 if ($hash{UTSWBlacklist});
    my $exacaf = '';
    if ($hash{AC_POPMAX} && $hash{AN_POPMAX}) {
      @exacs = split(/,/,$hash{AC_POPMAX});
      my $ac = 0;
      foreach $val (@exacs) {
	$ac += $val if ($val =~ m/^\d+$/);
      }
      @exans = split(/,/,$hash{AN_POPMAX});
      my $an = 0;
      foreach $val (@exans) {
	$an += $val if ($val =~ m/^\d+$/);
      }
      $exacaf = sprintf("%.4f",$ac/$an) if ($ac > 0 && $an > 10);
    }
    unless ($exacaf eq '' || $exacaf <= 0.01) {
      $fail{'COMMON'} = 1;
    }
    my $cosmicsubj = 0;
    if ($hash{CNT}) {
      my @cosmicct = split(/,/,$hash{CNT}); 
      foreach $val (@cosmicct) {
	$cosmicsubj += $val if ($val =~ m/^\d+$/);
      }
    }
    my @maf;
    my @dp;
    my @ao;
    my @ad;
    my @genotypes = @gts;
    my @deschead = split(/:/,$format);
  F1:foreach my $subjid (@gtheader) {
      my $allele_info = shift @gts;
      @ainfo = split(/:/, $allele_info);
      my %gtinfo = ();
      my @mutallfreq = ();
      foreach my $k (0..$#ainfo) {
	$gtinfo{$deschead[$k]} = $ainfo[$k];
	$hash{$deschead[$k]} = $ainfo[$k] if ($subjid eq $tumorid);
	$hash{'Normal'.$deschead[$k]} = $ainfo[$k] if ($subjid ne $tumorid);
      }
      my @altct = split(/,/,$gtinfo{AO});
      if ($gtinfo{DP} < 10) {
	if ($subjid eq $tumorid) {
	  next W1;
	}else {
	  delete $somline{$chrom}{$pos};
	  next W1;
	}
      }
      foreach  my $act (@altct) {
	push @mutallfreq, sprintf("%.4f",$act/$gtinfo{DP});
      }
      push @dp, $gtinfo{DP};
      push @ad, $gtinfo{AD};
      push @maf, \@mutallfreq;
      my @sortao = sort {$b <=> $a} @altct;
      push @ao, $sortao[0];
    }
    if ($gtheader[1] eq $tumorid) {
      @maf = reverse(@maf);
      @dp = reverse(@dp);
      @ad = reverse(@ad);
      @ao = reverse(@ao);
      @genotypes = reverse(@genotypes);
    }
    next if ($dp[0] < 20);
    $hash{AF} = join(",",@{$maf[0]});
    $hash{NormalAF} =  join(",",@{$maf[1]});
    my $newgt = $genotypes[0];
    $fail{'LowDepth'} = 1 if ($dp[0] < 20);
    my @callers = split(/,/,$hash{CallSet});
    if ($id =~ m/COS/ && $cosmicsubj >= 5) {
      $fail{'LowAltCt'} = 1 if ($ao[0] < 3);
      $fail{'LowMAF'} = 1 if ($maf[0][0] < 0.01);
    }else {
      $fail{'OneCaller'} = 1 if (scalar(@callers) < 2);
      $fail{'LowAltCt'} = 1 if ($ao[0] < 8);
      $fail{'LowMAF'} = 1 if ($maf[0][0] < 0.05);
    }
    if ($rnaval{$chrom}{$pos}) {
      my ($rnaad,$rnamaf,$rnadp) = @{$rnaval{$chrom}{$pos}};
      $hash{RnaSeqValidation} = 1 if ($rnamaf && $rnamaf =~ m/\d+/ && $rnamaf > 0.01);
      $hash{RnaSeqAF} = $rnamaf;
      $hash{RnaSeqDP} = $rnadp;
      $hash{RnaSeqAD} = $rnaad;
      $hash{fpkm} = $rnadp;
    }
    delete $hash{SOMATIC};
    if ($maf[1][0] >= 0.30) {
      $hash{Germline} = 1;
    }elsif ($maf[0][0] < 0.05 && ($maf[1][0] > 0.01 || $maf[1][0]*5 > $maf[0][0])) {
      next;
    }elsif ($maf[1][0] >= 0.05 || $maf[1][0]*5 > $maf[0][0]) {
      #$hash{Somatic} = 0;
      $hash{'LowFreqNormalAF'} = 1;
    }else {
      $hash{Somatic} = 1;
      $hash{SomaticCallSet}=$hash{CallSet};
    }
    next unless ($hash{ANN});
    foreach $trx (split(/,/,$hash{ANN})) {
      my ($allele,$effect,$impact,$gene,$geneid,$feature,
	  $featureid,$biotype,$rank,$codon,$aa,$pos_dna,$len_cdna,
	  $cds_pos,$cds_len,$aapos,$aalen,$distance,$err) = split(/\|/,$trx);
      next unless ($impact =~ m/HIGH|MODERATE/);
      next unless $keep{$gene};
      push @aa, $aa if ($aa ne '');
      $keepforvcf = $gene;
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
    print OUT join("\t",$chrom, $pos,$id,$ref,$alt,$score,
		   $filter,$newannot,$format,$newgt),"\n";
    next W1;
  }
}
close IN;

open IN, "gunzip -c $inputdir\/$tumorid/$tumorid\.annot.vcf.gz|" or die $!;
my %done;
W1:while (my $line = <IN>) {
  chomp($line);
  if ($line =~ m/^#CHROM/) {
    my @header = split(/\t/,$line);
    ($chrom, $pos,$id,$ref,$alt,$score,
     $filter,$info,$format,@gtheader) = split(/\t/, $line);
  }
  if ($line =~ m/^#/) {
    if ($somatic eq 'no_normal') {
      print OUT $line,"\n";
    }
    next;
  }
  my ($chrom, $pos,$id,$ref,$alt,$score,
      $filter,$annot,$format,@gts) = split(/\t/, $line);
  next if ($ref =~ m/\./ || $alt =~ m/\./ || $alt=~ m/,X/);
  next if ($somline{$chrom}{$pos});
  my %hash = ();
  foreach $a (split(/;/,$annot)) {
    my ($key,$val) = split(/=/,$a);
    $hash{$key} = $val unless ($hash{$key});
  }
  if ($rnaval{$chrom}{$pos}) {
    my ($rnaad,$rnamaf,$rnadp) = @{$rnaval{$chrom}{$pos}};
    $hash{RnaSeqValidation} = 1 if ($rnamaf && $rnamaf =~ m/\d+/ && $rnamaf > 0.01);
    $hash{RnaSeqAF} = $rnamaf;
    $hash{RnaSeqDP} = $rnadp;
    $hash{RnaSeqAD} = $rnaad;
    $hash{fpkm} = $rnadp;
  }
  my %fail;
  $fail{'UTSWBlacklist'} = 1 if ($hash{UTSWBlacklist});
  my $exacaf = '';
  if ($hash{AC_POPMAX} && $hash{AN_POPMAX}) {
    @exacs = split(/,/,$hash{AC_POPMAX});
    my $ac = 0;
    foreach $val (@exacs) {
      $ac += $val if ($val =~ m/^\d+$/);
    }
    @exans = split(/,/,$hash{AN_POPMAX});
    my $an = 0;
    foreach $val (@exans) {
      $an += $val if ($val =~ m/^\d+$/);
    }
    $exacaf = sprintf("%.4f",$ac/$an) if ($ac > 0 && $an > 10);
  }
  unless ($exacaf eq '' || $exacaf <= 0.01) {
    $fail{'COMMON'} = 1;
  }
  my @deschead = split(/:/,$format);
  my $allele_info = $gts[0];
  @ainfo = split(/:/, $allele_info);
  foreach my $k (0..$#ainfo) {
    $hash{$deschead[$k]} = $ainfo[$k];
  }
  my $cosmicsubj = 0;
  if ($hash{CNT}) {
    my @cosmicct = split(/,/,$hash{CNT}); 
    foreach $val (@cosmicct) {
      $cosmicsubj += $val if ($val =~ m/^\d+$/);
    }
  }
  my ($refct,@acts) = split(/,/,$hash{AD});
  my $totalaltct = 0;
  my @altnts = split(/,/,$alt);
  my @newalts;
  my @altct;
  my @normaltct;
  my @normmaf;
  my %normalcts;
  if ($normals{$chrom}{$pos}) {
    %normalcts = %{$normals{$chrom}{$pos}{BaseCts}};
    $hash{NormalDP} = $normals{$chrom}{$pos}{NormalDP};
    my $refct = 0;
    my $refdepth = 0;
    my $refaf = 0;
    ($refct,$refdepth,$refaf) = @{$normalcts{$ref}} if ($normalcts{$ref});
    push @normaltct, $refct;
  }
  foreach  my $i (0..$#acts) {
    if ($acts[$i] > 2) {
      if ($hash{NormalDP}) {
	my $normbasect = 0;
	my $normmaf = 0;
	($normbasect,$normdepth,$normmaf) = @{$normalcts{$altnts[$i]}} if ($normalcts{$altnts[$i]}) ;
	push @normaltct,$normbasect;
	push @normmaf,  $normmaf;
      }
      push @newalts, $altnts[$i];
      push @altct, $acts[$i];
      $totalaltct += $acts[$i];
    }
  }
  $hash{NormalAF} = join(",",@normmaf) if scalar(@normmaf> 0);
  $hash{NormalAD} = join(",",@normaltct) if scalar(@normmaf> 0);
  next unless ($altct[0] && $altct[0] > 2);
  if ($hash{DP} =~ m/,/) {
    $hash{DP} = $totalaltct+$hash{RO};
  }
  $hash{AD} = join(",",$hash{RO},@altct);
  next unless ($hash{DP});
  my @mutallfreq;
  foreach  my $act (@altct) {
    push @mutallfreq, sprintf("%.4f",$act/$hash{DP});
  }
  my @sortao = sort {$b <=> $a} @altct;
  $hash{AF} = join(",",@mutallfreq);
  if ($normals{$chrom}{$pos}) {
    if ($normmaf[0] >= 0.30) {
      $hash{Germline} = 1;
    }elsif ($mutallfreq[0] < 0.05 && ($normmaf[0] > 0.01 || $normmaf[0]*5 > $mutallfreq[0])) {
      next;
    }elsif ($normmaf[0] >= 0.05 || $normmaf[0]*5 > $mutallfreq[0]) {
      #$hash{Somatic} = 0;
      $hash{'LowFreqNormalAF'} = 1;
    }
  }
  if ($hash{DP} < 20) {
    next;
  }
  my @callers = split(/,/,$hash{CallSet});
  if ((grep(/hotspot/,@callers) || $id =~ m/COS/) && $cosmicsubj >= 5) {
    $fail{'LowAltCt'} = 1 if ($altct[0] < 3);
    $fail{'LowMAF'} = 1 if ($mutallfreq[0] < 0.01);
  }else {
    $fail{'OneCaller'} = 1 if (scalar(@callers) < 2);
    $fail{'LowAltCt'} = 1 if ($altct[0] < 8);
    $fail{'LowMAF'} = 1 if ($mutallfreq[0] < 0.05);
  }
  my $keepforvcf = 0;
  my @aa;
  next unless ($hash{ANN});
  foreach $trx (split(/,/,$hash{ANN})) {
    my ($allele,$effect,$impact,$gene,$geneid,$feature,
	$featureid,$biotype,$rank,$codon,$aa,$pos_dna,$len_cdna,
	$cds_pos,$cds_len,$aapos,$aalen,$distance,$err) = split(/\|/,$trx);
    next unless $keep{$gene};
    next unless ($impact =~ m/HIGH|MODERATE/);
    push @aa, $aa if ($aa ne '');
    $keepforvcf = $gene;
  }
  next unless $keepforvcf;
  my @fail = sort {$a cmp $b} keys %fail;
  if (scalar(@fail) < 1) {
    $filter = 'PASS';
  }elsif (scalar(@fail) > 0) {
    $filter = join(";", 'FailedQC',@fail);
  }else {
    next;
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
  print OUT join("\t",$chrom, $pos,$id,$ref,join(",",@newalts),$score,
		 $filter,$newannot,$format,$allele_info),"\n";
}

system("vcf-sort $inputdir\/$tumorid\.final.vcf |grep -v LowFreqNormalAF |bgzip  > $inputdir\/$subject\.utsw.vcf.gz");
system("vcf-sort $inputdir\/$tumorid\.final.vcf |grep '#\\|LowFreqNormalAF' |bgzip  > $inputdir\/$subject\.LowFreqNormal.vcf.gz");
system("vcf-sort $inputdir\/$tumorid\.final.vcf |grep -v 'FailedQC' |bgzip  > $inputdir\/$subject\.PASS.vcf.gz");

system("rm $inputdir\/$tumorid\.final.vcf");


sub log2 {
  my $n = shift;
  return log($n)/log(10);
}
