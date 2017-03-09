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
my $inputdir = '/project/PHG/PHG_Clinical/validation/$subject';
system("tabix $inputdir\/$tumorid/$tumorid\.annot.vcf.gz");
if ($rnaseqid ne 'no_rnaseq') {
  open RNACT, "<$inputdir\/$rnaseqid\/$rnaseqid\.ct" or die $!;
  while (my $line = <RNACT>) {
    chomp($line);
    next if ($line =~ m/^#|Geneid/);
    my ($geneid,$chr,$start,$end,$strand,$length,$ct) = split(/\t/,$line);
    $rnaseqct{$geneid} = $ct;
    $rnaseqct{'total'} += $ct;
  }
  close RNACT;
  open RNACT, "<$inputdir\/$rnaseqid\/$rnaseqid\.ct" or die $!;
  while (my $line = <RNACT>) {
    chomp($line);
    next if ($line =~ m/^#|Geneid/);
    my ($geneid,$chr,$start,$end,$strand,$length,$ct) = split(/\t/,$line);
    $rnaseqct{$geneid} = $ct;
    $rnaseqct{'total'} += $ct;
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
  open IN, "gunzip -c $inputdir\/somatic/$somatic\.annot.vcf.gz" or die $!;
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
    next unless ($exacaf eq '' || $exacaf <= 0.01);
    my @maf;
    my @dp;
    my @ao;
    my @newgt = @gt[0..1];
    my @deschead = split(/:/,$format);
  F1:foreach $subjid (@gtheader) {
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
    if ($gtheader[1] eq $tid) {
      @maf = reverse(@maf);
      @dp = reverse(@dp);
      @ao = reverse(@ao);
      @newgt = @gt[2..3];
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
    my $keepforvcf = 0;
    foreach $trx (split(/,/,$hash{ANN})) {
      my ($allele,$effect,$impact,$gene,$geneid,$feature,
	  $featureid,$biotype,$rank,$codon,$aa,$pos_dna,$len_cdna,
	  $cds_pos,$cds_len,$aapos,$aalen,$distance,$err) = split(/\|/,$trx);
      $aa =~ s/p\.//;
      $aa =~ s/($check)/$aa_hash{$1}/g;;
      #next unless ($impact =~ m/HIGH|MODERATE/ && $aa ne '');
      next unless $keep{$gene};
      $aa =~ m/([A-Z]+)\d+([A-Z]+)/;
      $keepforvcf = 1;
    }
    if ($keepforvcf) {
      $somline{$chrom}{$pos} = [$chrom, $pos,$id,$ref,$alt,$score,$filter,$annot,$format,@newgt];
      $somval{$chrom}{$pos} = $hash{CallSet};
    }
  }
}
close IN;

open IN, "gunzip -c $inputdir\/$tumorid/$tumorid\.annot.vcf.gz|" or die $!;
open OUT, ">$inputdir\/$tumorid\.final.vcf" or die $!;

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
  next unless ($exacaf eq '' || $exacaf <= 0.01);
  my @deschead = split(/:/,$format);
  my $allele_info = shift @gts;
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
  if (($hash{CallSet} =~ m/hotspot/ || $id =~ m/COS/) & $cosmicsubj >= 5) {
    next if ($dp < 20 || $altct < 3);
    next if ($maf < 0.01);
  }else {
      next unless ($hash{CallSet} =~ m/,/);
      next if ($dp < 20 || $altct < 8);
      next if ($maf < 0.05);
  }
  my $keepforvcf = 0;
  foreach $trx (split(/,/,$hash{ANN})) {
    my ($allele,$effect,$impact,$gene,$geneid,$feature,
	$featureid,$biotype,$rank,$codon,$aa,$pos_dna,$len_cdna,
	$cds_pos,$cds_len,$aapos,$aalen,$distance,$err) = split(/\|/,$trx);
    $aa =~ s/p\.//;
    $aa =~ s/($check)/$aa_hash{$1}/g;;
    next unless ($impact =~ m/HIGH|MODERATE/ && $aa ne '');
    next unless $keep{$gene};
    $done{$chrom}{$pos} = 1;
    $keepforvcf = $gene;
  }
  if ($keepforvcf) {
    $filter = 'PASS';
    if ($rnaseqct{$keepforvcf} > 10) {
      $annot .= ';logcpm='.log2(1000000*$rnaseqct{$keepforvcf}/$rnaseqct{'total'});
    } if ($logcpm{$keepforvcf}) {
      $annot .= ';fpkm='.$fpkm{$keepforvcf};
    } if ($rnaval{$chrom}{$pos}) {
      $annot .= ';RnaSeqValidation';
    } if ($somval{$chrom}{$pos}) {
      $annot .= ';Somatic;SomaticCallSet='.$somval{$chrom}{$pos};
    }
    print OUT join("\t",$chrom, $pos,$id,$ref,$alt,$score,$filter,$annot,$format,@gts),"\n";
  }
}
close OUT;

sub log2 {
  my $n = shift;
  return log($n)/log(10);
}
