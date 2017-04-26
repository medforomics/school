#!/usr/bin/perl -w
#integrate_datasets.pl

open OM, "</project/shared/bicf_workflow_ref/GRCh38/oncomine.genelist.txt" or die $!;
while (my $line = <OM>) {
    chomp($line);
    $oncomine{$line} = 'oncomine';
  }
close OM;
open OM, "</project/shared/bicf_workflow_ref/GRCh38/cosmic_consenus.genelist.txt" or die $!;
while (my $line = <OM>) {
    chomp($line);
    $cosmic{$line} = 'cosmic_consensus';
  }
close OM;
open OM, "</project/shared/bicf_workflow_ref/GRCh38/panel1385.genelist" or die $!;
while (my $line = <OM>) {
    chomp($line);
    $keep{$line} = 1;
}
close OM;

my ($subject,$samplename,$tumorid,$somatic,$rnaseqid) = @ARGV;
my %rnaseqct;
my $inputdir = "/project/PHG/PHG_Clinical/validation/$subject";
system("tabix -f $inputdir\/$tumorid/$tumorid\.annot.vcf.gz");
if ($rnaseqid ne 'no_rnaseq') {
  open RNACT, "<$inputdir\/$rnaseqid\/$rnaseqid\.cts" or die $!;
  while (my $line = <RNACT>) {
    chomp($line);
    next if ($line =~ m/^#|Geneid/);
    my ($geneid,$chr,$start,$end,$strand,$length,$ct) = split(/\t/,$line);
    $rnaseqct{$geneid} = $ct if ($ct);
    $rnaseqct{'total'} += $ct if ($ct);
  }
  close RNACT;
  open RNACT, "<$inputdir\/$rnaseqid\/$rnaseqid\.fpkm.txt" or die $!;
  while (my $line = <RNACT>) {
    chomp($line);
    next if ($line =~ m/^#|Geneid/);
    my ($ensid,$gene,$chr,$strand,$start,$end,$cov,$fpkm,$tmp) = split(/\t/,$line);
    $fpkm{$gene} = $fpkm if ($fpkm > 1);
  }
  close RNACT;
  
  system("zcat $inputdir\/$rnaseqid\/$rnaseqid\.annot.vcf.gz | perl -p -e 's/^/chr/' |perl -p -e 's/chr#/#/' |bgzip > $inputdir\/$rnaseqid\/$rnaseqid\.annot.chr.vcf.gz");
  system("tabix $inputdir\/$rnaseqid\/$rnaseqid\.annot.chr.vcf.gz");
  system("bcftools isec -p $inputdir\/rnaoverlap --nfiles =2 $inputdir\/$rnaseqid\/$rnaseqid\.annot.chr.vcf.gz $inputdir\/$tumorid/$tumorid\.annot.vcf.gz");
  open RNAOV, "<$inputdir\/rnaoverlap/sites.txt" or die $!;
  while (my $line = <RNAOV>) {
    chomp($line);
    my ($chr,$pos,$ref,$alt,$overlap) = split(/\t/,$line);
    $rnaval{$chr}{$pos} = 1;
  }
}
if ($somatic ne 'no_normal') {
  open IN, "gunzip -c $inputdir\/somatic/$somatic\.annot.vcf.gz |" or die $!;
 W1:while (my $line = <IN>) {
    chomp($line);
    if ($line =~ m/^#CHROM/) {
      my @header = split(/\t/,$line);
      ($chrom, $pos,$id,$ref,$alt,$score,
       $filter,$info,$format,@gtheader) = split(/\t/, $line);
    }
    if ($line =~ m/^#/) {
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
    #next unless ($exacaf eq '' || $exacaf <= 0.01);
    my @maf;
    my @dp;
    my @ao;
    my $newgt = $gts[0];
    if ($gtheader[1] eq $tumorid) {
      $newgt = $gts[1];
    }
    my @deschead = split(/:/,$format);
  F1:foreach my $subjid (@gtheader) {
      my $allele_info = shift @gts;
      @ainfo = split(/:/, $allele_info);
      my %gtinfo = ();
      foreach my $k (0..$#ainfo) {
	$gtinfo{$deschead[$k]} = $ainfo[$k];
      }
      my $altct = 0;
      foreach  my $act (split(/,/,$gtinfo{AO})) {
	$altct += $act;
      }
      $gtinfo{AO} = $altct;
      next W1 if ($gtinfo{DP} < 10);
      push @dp, $gtinfo{DP};
      push @maf, sprintf("%.4f",$gtinfo{AO}/$gtinfo{DP});
      push @ao, $gtinfo{AO};
    }
    if ($gtheader[1] eq $tumorid) {
      @maf = reverse(@maf);
      @dp = reverse(@dp);
      @ao = reverse(@ao);
    }
    my $cosmicsubj = 0;
    @cosmicct = split(/,/,$hash{CNT}) if $hash{CNT};
    foreach $val (@exacs) {
      $cosmicsubj += $val if ($val =~ m/^\d+$/);
    }
    next if ($maf[1] > 0.005 && $maf[1]*5 > $maf[0]);
    next if ($maf[1] > 0.005);
    next if ($maf[0] < 0.01);
    if ($id =~ m/COSM/ && $cosmicsubj >=5) {
      next if ($ao[0] < 3);
      next if ($maf[0] < 0.01);
    }else {
      next if ($ao[0] < 8);
      next if ($hash{CallSet} !~ m/,/);
      next if ($maf[0] < 0.05);
    }
    if ($rnaval{$chrom}{$pos}) {
      $annot .= ';RnaSeqValidation';
    }
    $annot .= ';Somatic;SomaticCallSet='.$hash{CallSet};
    foreach $trx (split(/,/,$hash{ANN})) {
      my ($allele,$effect,$impact,$gene,$geneid,$feature,
	  $featureid,$biotype,$rank,$codon,$aa,$pos_dna,$len_cdna,
	  $cds_pos,$cds_len,$aapos,$aalen,$distance,$err) = split(/\|/,$trx);
      next unless $keep{$gene};
      if ($rnaseqct{$gene}) {
	$annot .= ';logcpm='.sprintf("%.1f",log2(1000000*$rnaseqct{$gene}/$rnaseqct{'total'}));
      } if ($fpkm{$gene}) {
	$annot .= ';fpkm='.sprintf("%.1f",$fpkm{$gene});
      } if ($oncomine{$gene} || $cosmic{$gene}) {
	$annot .= ';CancerGene';
      }
      if (($exacaf eq '' || $exacaf <= 0.01) && $impact =~ m/HIGH|MODERATE/) {
	$somline{$chrom}{$pos} = [$chrom, $pos,$id,$ref,$alt,$score,
				  'PASS',$annot,$format,$newgt];
      }
      $somval{$chrom}{$pos} = $hash{CallSet};
      next W1;
    }
  }
}
close IN;

open IN, "gunzip -c $inputdir\/$tumorid/$tumorid\.annot.vcf.gz|" or die $!;
open OUT, ">$inputdir\/$tumorid\.final.vcf" or die $!;
my %done;
while (my $line = <IN>) {
  chomp($line);
  if ($line =~ m/^#CHROM/) {
    my @header = split(/\t/,$line);
    ($chrom, $pos,$id,$ref,$alt,$score,
     $filter,$info,$format,@gtheader) = split(/\t/, $line);
  }
  if ($line =~ m/^#/) {
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
  my %fail;
  unless ($exacaf eq '' || $exacaf <= 0.01) {
    $fail{'COMMON'} = 1;
  }
  my @deschead = split(/:/,$format);
  my $allele_info = $gts[0];
  @ainfo = split(/:/, $allele_info);
  my %gtinfo = ();
  foreach my $k (0..$#ainfo) {
      $gtinfo{$deschead[$k]} = $ainfo[$k];
  }
  my $altct = (split(/,/,$gtinfo{AO}))[0];
  foreach  my $act (split(/,/,$gtinfo{AO})) {
    $altct += $act;
  }
  my $dp = $gtinfo{DP};
  if ($gtinfo{DP} =~ m/,/) {
    $dp = $altct+$gtinfo{RO};
  }
  $maf = sprintf("%.4f",$altct/$dp);
  my $cosmicsubj = 0;
  @cosmicct = split(/,/,$hash{CNT}) if $hash{CNT};
  foreach $val (@exacs) {
    $cosmicsubj += $val if ($val =~ m/^\d+$/);
  }
  if ($dp < 20) {
    $fail{'LowDepth'} = 1;
  }
  if (($hash{CallSet} =~ m/hotspot/ || $id =~ m/COS/) && $cosmicsubj >= 5) {
    $fail{'LowAltCt'} = 1 if ($altct < 3);
    $fail{'LowMAF'} = 1 if ($maf < 0.01);
  }else {
    $fail{'OneCaller'} = 1 if ($hash{CallSet} =~ m/,/);
    $fail{'LowAltCt'} = 1 if ($dp < 20 || $altct < 8);
    $fail{'LowMAF'} = 1 if ($maf < 0.05);
  }
  my $keepforvcf = 0;
  my @aa;
  foreach $trx (split(/,/,$hash{ANN})) {
    my ($allele,$effect,$impact,$gene,$geneid,$feature,
	$featureid,$biotype,$rank,$codon,$aa,$pos_dna,$len_cdna,
	$cds_pos,$cds_len,$aapos,$aalen,$distance,$err) = split(/\|/,$trx);
    next unless ($impact =~ m/HIGH|MODERATE/);
    next unless $keep{$gene};
    push @aa, $aa if ($aa ne '');
    $keepforvcf = $gene;
  }
  $fail{'NoAAChange'} = 1 if (scalar(@aa) < 1);
  my @fail = keys %fail;
  if ($keepforvcf && scalar(@fail) < 1) {
    $filter = 'PASS';
  }elsif (scalar(@fail) > 0) {
    $filter = join(";", 'REMOVE',@fail);
  }else {
    next;
  }
  $done{$chrom}{$pos} = 1;
  if ($keepforvcf) {
    if ($rnaseqct{$keepforvcf} > 10) {
      $annot .= ';logcpm='.sprintf("%.1f",log2(1000000*$rnaseqct{$keepforvcf}/$rnaseqct{'total'}));
    } if ($fpkm{$keepforvcf}) {
      $annot .= ';fpkm='.sprintf("%.1f",$fpkm{$keepforvcf});
    } if ($oncomine{$keepforvcf} || $cosmic{$keepforvcf}) {
      $annot .= ';CancerGene';
    }
  } if ($rnaval{$chrom}{$pos}) {
    $annot .= ';RnaSeqValidation';
  } if ($somval{$chrom}{$pos}) {
    $annot .= ';Somatic;SomaticCallSet='.$somval{$chrom}{$pos};
  }
  print OUT join("\t",$chrom, $pos,$id,$ref,$alt,$score,$filter,$annot,$format,@gts),"\n";
}

foreach my $chrom (keys %somline) {
  foreach my $pos (keys %{$somline{$chrom}}) {
    next if $done{$chrom}{$pos};
    print OUT join("\t",@{$somline{$chrom}{$pos}}),"\n";
  }
}
close OUT;

system("vcf-sort $inputdir\/$tumorid\.final.vcf | bgzip > $inputdir\/$tumorid\.final.vcf.gz");
system("rm $inputdir\/$tumorid\.final.vcf");
system("rm -fr $inputdir\/rnaoverlap");

sub log2 {
  my $n = shift;
  return log($n)/log(10);
}
