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
open OM, "</project/shared/bicf_workflow_ref/GRCh38/panel1385.genelist.txt" or die $!;
while (my $line = <OM>) {
    chomp($line);
    $keep{$line} = 1;
}
close OM;

my @normals = @ARGV;
foreach my $norm (@normals) {
  chomp($norm);
  open IN, "gunzip -c $norm |" or die $!;
  $outfile = $norm;
  $outfile =~ s/annot.vcf.gz/filt.vcf/;
  open OUT, ">$outfile" or die $!;
 W1:while (my $line = <IN>) {
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
    my @cts = split(/,/,$gtinfo{AD});
    my $refct = shift @cts;
    my $altct = 0;
    foreach  my $act (@cts) {
      $altct += $act;
    }
    my $dp = $altct+$refct;
    $dp = $gtinfo{DP} unless ($dp);
    next W1 unless ($dp && $dp > 9);
    $maf = sprintf("%.4f",$altct/$dp);
    my $cosmicsubj = 0;
    @cosmicct = split(/,/,$hash{CNT}) if $hash{CNT};
    foreach $val (@cosmicct) {
      $cosmicsubj += $val if ($val =~ m/^\d+$/);
    }
    if ($id =~ m/COSM/ && $cosmicsubj >=5) {
      next if ($altct < 3 || $maf < 0.01);
    }else {
      next if ($altct < 8 || $maf < 0.05);
      next if ($hash{CallSet} !~ m/,/);
    }
    my $keep = 0;
    foreach $trx (split(/,/,$hash{ANN})) {
      my ($allele,$effect,$impact,$gene,$geneid,$feature,
	  $featureid,$biotype,$rank,$codon,$aa,$pos_dna,$len_cdna,
	  $cds_pos,$cds_len,$aapos,$aalen,$distance,$err) = split(/\|/,$trx);
      next unless $keep{$gene};
      if ($impact =~ m/HIGH|MODERATE/) {
	$keep = 1;
      }
    }
    print OUT $line,"\n" if ($keep);
  }
}
close IN;
close OUT;
