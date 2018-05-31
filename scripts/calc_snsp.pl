#!/usr/bin/perl -w
#calc_snp.pl

my %tp;
my %fp;
my %fn;
my $input = shift @ARGV;
my $fn = $input;
$fn =~ s/eval.vcf/utswcoding.multiinter.bed.gz/;
my $run = $input;
$run =~ s/\.eval.vcf//;
my @allcallers = ('strelka2','sam','ssvar','platypus','gatk');
foreach ((@allcallers,'union','genomeseer')) {
  $fn{'snp'}{$_} = 0;
  $fn{'indel'}{$_} = 0;
}
open FP, ">fp.vcf" or die $!;
open FN, "gunzip -c $fn |" or die $!;
while (my $line = <FN>) {
  chomp($line);
  my ($chr,$pos,$end,$calls) = split(/\t/,$line);
  next if ($calls =~ m/union/);
  next unless ($calls =~ m/platinum/ && $calls =~ m/giab/);
  $vartype = 'snp';
  $size = $end-$pos;
  if ($size > 1) {
    $vartype = 'indel';
  }
  next if ($chr eq 'chr6' && $pos >= 28510120 && $pos <= 33480577);
  $fn{$vartype}{'union'} ++;
  $fn{$vartype}{'genomeseer'} ++;
  print $line,"\n";
  foreach (@allcallers) {
    $fn{$vartype}{$_} ++;
  }
}

open IN, "<$input" or die $!;
while (my $line = <IN>) {
  chomp($line);
  if ($line =~ m/^#/) {
    print FP $line,"\n";
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
  my $cosmicsubj = 0;
  @cosmicct = split(/,/,$hash{CNT}) if $hash{CNT};
  foreach $val (@cosmicct) {
    $cosmicsubj += $val if ($val =~ m/^\d+$/);
  }
  
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
  my $maf = sprintf("%.4f",$altct/$dp);
  my $vartype = 'snp';
  if ($hash{TYPE} ne 'snp') {
    next if $maf < 0.1;
    $vartype = 'indel';
  }
  #next if (exists $hash{FS} && $hash{FS} > 60);
  $is_gs = 1;
  $is_gs = 0 if ($hash{CallSet} !~ m/\|/ && $hash{CallSet} =~ m/hotspot/ && $maf > 0.1);
  if (($id =~ m/COS/) && $cosmicsubj >= 5) {
    $is_gs = 0 if ($dp < 20 || $altct < 3);
    $is_gs = 0 if ($maf < 0.01 && $vartype eq 'snp');
    $is_gs = 0 if ($maf < 0.03 && $vartype eq 'indel');
  }else {
    $is_gs = 0 unless ($hash{CallSet} =~ m/\|/);
    $is_gs = 0 if ($dp < 20 || $altct < 8);
    $is_gs = 0 if ($maf < 0.05 && $vartype eq 'snp');
    $is_gs = 0 if ($maf < 0.1 && $vartype eq 'indel');
  }
  %chash = map {(split(/\//))[0]=>1} split(/\|/,$hash{CallSet});
  if ($hash{PlatRef} =~ m/giab|platinum/) {
    if ($is_gs){
      $tp{$vartype}{'genomeseer'} ++;
    }elsif($is_gs && $hash{PlatRef} =~ m/giab/ && $hash{PlatRef} =~ m/platinum/) {
      $fn{$vartype}{'genomeseer'} ++;
    }
    $tp{$vartype}{'union'} ++;
    foreach (@allcallers) {
      if ($chash{$_}){
	$tp{$vartype}{$_} ++;
      }elsif($hash{PlatRef} =~ m/giab/ && $hash{PlatRef} =~ m/platinum/) {
	$fn{$vartype}{$_} ++;
      }
    }
  }else {
    if ($is_gs){
      $fp{$vartype}{'genomeseer'} ++;
      print FP $line,"\n";
    }
    $fp{$vartype}{'union'} ++;
    foreach $caller (keys %chash) {
      $fp{$vartype}{$caller} ++;
    }
  }
}

foreach ((@allcallers,'genomeseer','union')) {
  $fp{'snp'}{$_} = 0 unless $fp{'snp'}{$_};
  $fp{'indel'}{$_} = 0 unless $fp{'indel'}{$_};
  $fn{'indel'}{$_} = 0 unless $fn{'indel'}{$_};
  next unless ($tp{'snp'}{$_});
  $sn_snp = sprintf("%.1f",100*$tp{'snp'}{$_}/($tp{'snp'}{$_}+$fn{'snp'}{$_}));
  $sn_indel = sprintf("%.1f",100*$tp{'indel'}{$_}/($tp{'indel'}{$_}+$fn{'indel'}{$_}));
  $sp = sprintf("%.1f",100*($tp{'snp'}{$_}+$tp{'indel'}{$_})/($tp{'snp'}{$_}+$fp{'snp'}{$_}+$tp{'indel'}{$_}+$fp{'indel'}{$_}));
  print join("\t",$run,$_,$tp{'snp'}{$_},$fp{'snp'}{$_},$fn{'snp'}{$_},
	     $tp{'indel'}{$_},$fp{'indel'}{$_},$fn{'indel'}{$_},
	     $sn_snp,$sn_indel,$sp),"\n";
}
