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
my @allcallers = ('gatk','sam','ssvar','platypus');
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
  foreach (@allcallers) {
    $fn{$vartype}{$_} ++;
  }
}
open IN, "<$input" or die $!;
while (my $line = <IN>) {
  chomp($line);
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
  if (exists $hash{TYPE} and $hash{TYPE} ne 'snp') {
    next if $maf < 0.1;
    $vartype = 'indel';
  }
  $is_gs = 1;
  if (($hash{CallSet} =~ m/hotspot/ || $id =~ m/COS/) && $cosmicsubj >= 5) {
    $is_gs = 0 if ($dp < 25 || $altct < 3);
    $is_gs = 0 if ($maf < 0.01 && $vartype eq 'snp');
    $is_gs = 0 if ($maf < 0.03 && $vartype eq 'indel');
  }else {
    $is_gs = 0 unless ($hash{CallSet} =~ m/sam/ || $hash{CallSet} =~ m/,/);
    $is_gs = 0 if ($dp < 25 || $altct < 8);
    $is_gs = 0 if ($maf < 0.05 && $vartype eq 'snp');
    $is_gs = 0 if ($maf < 0.1 && $vartype eq 'indel');
  }
  %chash = map {$_=>1} split(/,/,$hash{CallSet});
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
      print $line,"\n";
    }
    $fp{$vartype}{'union'} ++;
    foreach $caller (keys %chash) {
      $fp{$vartype}{$caller} ++;
    }
  }
}

open OUT, ">$run\.calcsnps.summary" or die $!;
print OUT join ("\t","Sample_Name","Caller","TP_SNP","FP_SNP","FN_SNP","SNP_Sensitivity",
                "SNP_PPV","TP_INDEL","FP_INDEL","FN_INDEL","INDEL_Sensitivity","INDEL_PPV",
                "TP_Variant,FP_Variant,FN_Variant,Variant_Sensitivity,Variant_PPV"),"\n";
foreach ((@allcallers,'genomeseer','union')) {
  $fp{'snp'}{$_} = 0 unless $fp{'snp'}{$_};
  $fp{'indel'}{$_} = 0 unless $fp{'indel'}{$_};
  $fn{'indel'}{$_} = 0 unless $fn{'indel'}{$_};
  my $merged_tp = $tp{'snp'}{$_}+$tp{'indel'}{$_};
  my $merged_fp = $fp{'snp'}{$_}+$fp{'indel'}{$_};
  my $merged_fn = $fn{'snp'}{$_}+$fn{'indel'}{$_};
  print OUT join("\t",$run,$_,$tp{'snp'}{$_},$fp{'snp'}{$_},$fn{'snp'}{$_},
	     sprintf("%.1f",100*$tp{'snp'}{$_}/($tp{'snp'}{$_}+$fn{'snp'}{$_})),
	     sprintf("%.1f",100*$tp{'snp'}{$_}/($tp{'snp'}{$_}+$fp{'snp'}{$_})),
	     $tp{'indel'}{$_},$fp{'indel'}{$_},$fn{'indel'}{$_},
	     sprintf("%.1f",100*$tp{'indel'}{$_}/($tp{'indel'}{$_}+$fn{'indel'}{$_})),
	     sprintf("%.1f",100*$tp{'indel'}{$_}/($tp{'indel'}{$_}+$fp{'indel'}{$_})),
             $merged_tp,$merged_fp,$merged_fn, sprintf("%.1f",100*$merged_tp/($merged_tp+$merged_fn)),
             sprintf("%.1f",100*$merged_tp/($merged_tp+$merged_fp))),"\n";
}
