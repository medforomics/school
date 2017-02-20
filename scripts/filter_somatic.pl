#!/usr/bin/perl -w
#use strict;

my %aa_hash=(Ala=>'A',Arg=>'R',Asn=>'N',Asp=>'D',Cys=>'C',
	     Glu=>'E',Gln=>'Q',Gly=>'G',His=>'H',Ile=>'I',
	     Leu=>'L',Lys=>'K',Met=>'M',Phe=>'F',Pro=>'P',
	     Ser=>'S',Thr=>'T',Trp=>'W',Tyr=>'Y',Val=>'V');

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
my $check = join '|', keys %aa_hash;
my $prefix = shift @ARGV;
my $tid = shift @ARGV;
my $input = "$prefix\.annot.vcf.gz" or die $!;
open OUT, ">$prefix\.coding_variants.txt" or die $!;
open SOM, ">$prefix\.somatic.vcf" or die $!;

open IN, "gunzip -c $input|" or die $!;
while (my $line = <IN>) {
  chomp($line);
  if ($line =~ m/^#CHROM/) {
    my @header = split(/\t/,$line);
    ($chrom, $pos,$id,$ref,$alt,$score,
     $filter,$info,$format,@gtheader) = split(/\t/, $line);
  }
  if ($line =~ m/^#/) {
    print SOM $line,"\n";
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
  my @maf;
  my @dp;
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
    push @dp, $gtinfo{DP};
    push @maf, sprintf("%.4f",$gtinfo{AO}/$gtinfo{DP});
  }
  if ($gtheader[1] eq $tid) {
    @maf = reverse(@maf);
    @dp = reverse(@dp);
  }
  next if ($dp[0] < 100 || $dp[1] < 25);
  next if ($hash{CallSet} !~ m/,/);
  next if ($maf[1] > 0.005);
  next if ($maf[0] < 0.01);
  next if ($maf[1] > 0 && $maf[1]*10 > $maf[0]);
  my $exacaf = '';
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
  @cosmicct = split(/,/,$hash{CNT});
  my $cosmicsubj = 0;
  foreach $val (@exacs) {
    $cosmicsubj += $val if ($val =~ m/^\d+$/);
  }
  $exacaf = sprintf("%.4f",$ac/$an) if ($ac > 0 && $an > 10);
  my $keepforvcf = 0;
  next unless ($exacaf eq '' || $exacaf <= 0.05);
  foreach $trx (split(/,/,$hash{ANN})) {
    my ($allele,$effect,$impact,$gene,$geneid,$feature,
	$featureid,$biotype,$rank,$codon,$aa,$pos_dna,$len_cdna,
	$cds_pos,$cds_len,$aapos,$aalen,$distance,$err) = split(/\|/,$trx);
    $aa =~ s/p\.//;
    $aa =~ s/($check)/$aa_hash{$1}/g;;
    next unless ($impact =~ m/HIGH|MODERATE/ && $aa ne '');
    next unless $keep{$gene};
    next if ($done{$chrom}{$pos});
    $aa =~ m/([A-Z]+)\d+([A-Z]+)/;
    $done{$chrom}{$pos} = 1;
    $oncomine{$gene} = '' unless $oncomine{$gene};
    $cosmic{$gene} = '' unless $cosmic{$gene};
    print OUT join("\t",$chrom,$pos,$id,@maf,@dp,$exacaf,$hash{CallSet},$cosmicsubj,$cosmic{$gene},$oncomine{$gene},$gene,$aa,$effect,$impact),"\n";
    $keepforvcf = 1;
  }
  if ($keepforvcf) {
    print SOM $line,"\n";
  }
}
close SOM;
