#!/usr/bin/perl -w
#integrate_datasets.pl

#module load vcftools/0.1.14 samtools/1.6 bedtools/2.26.0 
my $refdir = '/project/shared/bicf_workflow_ref/GRCh38/';
my $execdir = '/project/PHG/PHG_Clinical/clinseq_workflows/';
open OM, "<$refdir\/panel1385.genelist.txt" or die $!;
while (my $line = <OM>) {
  chomp($line);
  $keep{$line} = 1;
}

open OM, "<$refdir\/cancer.genelist.txt" or die $!;
while (my $line = <OM>) {
  chomp($line);
  $cgenelist{$line} = 1;
}


my ($subject,$tumorid,$tumorvcf,$somaticvcf,$rnaseqvcf,$rnaseqbam) = @ARGV;
if ($somaticvcf ne 'no_normal') {
  system("zcat $somaticvcf > somatic.vcf");
  system("$execdir/scripts/vcf2bed.pl somatic.vcf |cut -f 1,2,3 > somatic.bed");
  system("bedtools intersect -header -v -b somatic.bed -a $tumorvcf > tumoronly.vcf");
  system("vcf-shuffle-cols -t somatic.vcf tumoronly.vcf |bgzip > tumor.vcf.gz");
  system("vcf-concat $somaticvcf tumor.vcf.gz |vcf-sort |bgzip > somatic_germline.vcf.gz");
} else {
  system("ln -s $tumorvcf somatic_germline.vcf.gz");
}
system("tabix somatic_germline.vcf.gz");

if ($rnaseqvcf ne 'no_rnaseq') {
  system("zcat  $rnaseqvcf |perl -p -e 's/^/chr/g' | perl -p -e 's/^chr#/#/g' |bgzip > rnaseq.vcf.gz");
  system("tabix rnaseq.vcf.gz");
  system("vcf-merge somatic_germline.vcf.gz rnaseq.vcf.gz |bgzip > allvariants.vcf.gz");
} else {
  system("ln -s somatic_germline.vcf.gz allvariants.vcf.gz");
}

my $rnaseqid;
if ($rnaseqvcf ne 'no_rnaseq') {
  system("zcat somatic_germline.vcf.gz > alltumor.vcf");
  system("/project/PHG/PHG_Clinical/clinseq_workflows/scripts/vcf2bed.pl alltumor.vcf |cut -f 1,2,3 > tumor.bed");
  system("cat tumor.bed |perl -p -e 's/chr//g' > tumor.nochr.bed");
  system("/project/shared/bicf_workflow_ref/seqprg/bam-readcount/bin/bam-readcount -w 0 -q 0 -b 25 -l tumor.nochr.bed -f $refdir/\/hisat_genome.fa $rnaseqbam > rnaseq.bamreadct.txt") unless ( -e "rnaseq.bamreadct.txt");
  $idline = `samtools view -H ORD543-27-3044_T_RNA_panelrnaseq.final.bam |grep '^\@RG' |grep 'ID:'`;
  $idline =~ m/ID:(\S+)/;
  $rnaseqid = $1;
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
	$hash{$base} = $basect;
      }else {
	$hash{$base} = $basect if ($basect);
      }
    }
    $rnaval{'chr'.$chr}{$pos} = [\%hash,$depth];
  }
}
$somaticvcf =~ m/$tumorid\_(.+).annot/g;
$normalid = $1;

open OUT, ">$tumorid\.final.vcf" or die $!;
open PASS, ">$tumorid\.pass.vcf" or die $!;

open IN, "gunzip -c allvariants.vcf.gz |" or die $!;
W1:while (my $line = <IN>) {
  chomp($line);
  if ($line =~ m/^#CHROM/) {
    my @header = split(/\t/,$line);
    ($chrom, $pos,$id,$ref,$alt,$score,
     $filter,$info,$format,@gtheader) = split(/\t/, $line);
    print OUT join("\t",$chrom,$pos,$id,$ref,$alt,$score,
		   $filter,$info,$format,@gtheader),"\n";
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
  if ($hash{FS} && $hash{FS} > 60) {
    $fail{'StrandBias'} = 1;
  }
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
    next F1 unless ($gtinfo{$subjid}{DP} && $gtinfo{$subjid}{DP} >= 5);
    my @altct = split(/,/,$gtinfo{$subjid}{AO});
    foreach  my $act (@altct) {
	next if ($act eq '.');
      push @mutallfreq, sprintf("%.4f",$act/$gtinfo{$subjid}{DP});
    }
    $gtinfo{$subjid}{MAF} = \@mutallfreq;
  }
  next unless ($gtinfo{$tumorid}{DP} && $gtinfo{$tumorid}{DP} >= 20);
  @tumormaf = @{$gtinfo{$tumorid}{MAF}};
  @tumoraltct = split(/,/,$gtinfo{$tumorid}{AO});
  next if ($tumoraltct[0] eq '.');
  $hash{AF} = join(",",@tumormaf);
  my @callers = split(/,/,$hash{CallSet});
  if ($id =~ m/COS/ && $cosmicsubj >= 5) {
    $fail{'LowAltCt'} = 1 if ($tumoraltct[0] < 3);
    $fail{'LowMAF'} = 1 if ($tumormaf[0] < 0.01);
    $fail{'LowMAF'} = 1 if ($tumormaf[0] < 0.03 && $hash{TYPE} ne 'snp');
  }else {
    $fail{'OneCaller'} = 1 if (scalar(@callers) < 2);
    $fail{'LowAltCt'} = 1 if ($tumoraltct[0] < 8);
    $fail{'LowMAF'} = 1 if ($tumormaf[0] < 0.05);
    $fail{'LowMAF'} = 1 if ($tumormaf[0] < 0.10 && $hash{TYPE} ne 'snp');
  }
  $hash{SS} = 5;
  delete $hash{SOMATIC};
  if ($gtinfo{$normalid}) {
    next unless ($gtinfo{$normalid}{MAF});
    @normalmaf = @{$gtinfo{$normalid}{MAF}};
    if ($normalmaf[0] >= 0.30) {
      $hash{SS} = 1;
    }elsif ($tumormaf[0] < 0.05 && ($normalmaf[0] > 0.01 || $normalmaf[0]*5 > $tumormaf[0])) {
      next;
    }elsif ($normalmaf[0] >= 0.05 || $normalmaf[0]*5 > $tumormaf[0]) {
      $hash{SS} = 5;
      $hash{'LowFreqNormalAF'} = 1;
    }else {
      $hash{SS} = 2;
    }
  }
  my $rna_gtinfo = join(",",'.','.','.','.','.');
  if ($rnaval{$chrom}{$pos}) {
    next if ($gtinfo{$rnaseqid}{DP});
    my ($rnahashref,$rnadp) = @{$rnaval{$chrom}{$pos}};
    my %rnantct = %{$rnahashref};
    $rnantct{$ref} = 0 unless ($rnantct{$ref});
    $gtinfo{$rnaseqid}{RO} = $rnantct{$ref};
    my @altcts;
    foreach $altnt (split(/,/,$alt)) {
      my $ct = $rnantct{$altnt};
      $ct = 0 unless ($ct);
      push @altcts, $ct;
    }
    $gtinfo{$rnaseqid}{AO} = join(",",@altcts);
    $gtinfo{$rnaseqid}{GT} = '.';
    $gtinfo{$rnaseqid}{DP} = $rnadp;
    $gtinfo{$rnaseqid}{AD} = join(",",$gtinfo{$rnaseqid}{RO},@altcts);
  }
  $hash{RnaSeqValidation} = 1 if ($gtinfo{$rnaseqid}{AO} && 
				  $gtinfo{$rnaseqid}{AO} =~ m/^\d+/ &&
				  (split(/,/,$gtinfo{$rnaseqid}{AO}))[0] > 2);
  my @newgt;
  foreach $sample (@gtheader) {
    my @gtdata;
    foreach $gt (@deschead) {
      $gtinfo{$sample}{$gt} = '.' unless (exists $gtinfo{$sample}{$gt});
      push @gtdata, $gtinfo{$sample}{$gt};
    }
    push @newgt, join(":",@gtdata);
  }
  next unless ($hash{ANN});
  my $cancergene = 0;
  foreach $trx (split(/,/,$hash{ANN})) {
    my ($allele,$effect,$impact,$gene,$geneid,$feature,
	$featureid,$biotype,$rank,$codon,$aa,$pos_dna,$len_cdna,
	$cds_pos,$cds_len,$aapos,$aalen,$distance,$err) = split(/\|/,$trx);
    next unless ($impact =~ m/HIGH|MODERATE/);
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
		 $format,@newgt),"\n" if ($filter eq 'PASS');
  print OUT join("\t",$chrom, $pos,$id,$ref,$alt,$score,$filter,$newannot,
		 $format,@newgt),"\n" if ($filter eq 'PASS' || $id =~ m/COS/ || $cancergene);
  next W1;
}
close IN;

system("vcf-sort $tumorid\.final.vcf |bedtools intersect -header -a stdin -b $refdir\/UTSWV2.bed  |uniq |bgzip > $subject\.utsw.vcf.gz");
system("bedtools intersect -header -a $tumorid\.pass.vcf -b $refdir\/UTSWV2.bed |bgzip > $subject\.PASS.vcf.gz");
system("rm $tumorid\.final.vcf");

open TMB, ">$subject\.tmb.txt" or die $!;
$num_mutations = `zgrep -c -v "SS=2" $subject\.PASS.vcf.gz`;
chomp($num_mutations);
print TMB join("\n","Class,TMB",join(",",'',sprintf("%.2f",$num_mutations/4.6))),"\n";

sub log2 {
    my $n = shift;
    return log($n)/log(10);
}
