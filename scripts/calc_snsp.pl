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
open FP, ">$run\.fp.vcf" or die $!;
open GS, ">$run\.genomeseer.vcf" or die $!;
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
my %done;
open IN, "<$input" or die $!;
while (my $line = <IN>) {
  chomp($line);
  if ($line =~ m/^#/) {
    print FP $line,"\n";
    print GS $line,"\n";
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
  $cosmicsubj = $cosmicct[0];
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
  if (length($alt) == length($ref) && length($ref) > 1) {
    $vartype = 'indel';
  }
  foreach $ant (split(/,/,$alt)) {
    $vartype = 'indel' if (length($ant) != length($ref));
  }
  if ($hash{REFREP} && $hash{REFREP} =~ m/^\d$/) {
      $numreps= (split(/,/,$hash{REFREP}))[0];
      next if $numreps > 1;
  }
  my @callers;
  if ($hash{CallSet} && $hash{CallSet} =~ m/\// || $hash{SomaticCallSet} && $hash{SomaticCallSet} =~ m/\//) {
    my @callinfo ;
    @callinfo = split(/|/, $hash{CallSet}) if ($hash{CallSet});
    if ($hash{SomaticCallSet}) {
      @callinfo = split(/|/, $hash{SomaticCallSet}) if ($hash{SomaticCallSet});
    }
    foreach $cinfo (@callinfo) {
      my ($caller, $alt, @samafinfo) = split(/\//,$cinfo);
      push @callers, $caller;
    }
    $hash{CallSet} = join(",",@callinfo);
  } elsif ($hash{CallSet} && $hash{CallSet} =~ m/\|/ || $hash{SomaticCallSet} && $hash{SomaticCallSet} =~ m/\|/) {
    my @callinfo ;
    @callinfo = split(/,/, $hash{CallSet}) if ($hash{CallSet});
    if ($hash{SomaticCallSet}) {
      @callinfo = split(/,/, $hash{SomaticCallSet}) if ($hash{SomaticCallSet});
    }
    foreach $cinfo (@callinfo) {
      my ($caller, $alt, @samafinfo) = split(/\|/,$cinfo);
      push @callers, $caller;
    }
    $hash{CallSet} = join(",",@callinfo);
  } else {
    if ($hash{CallSet}) {
      @callers = (@callers, split(/\,/, $hash{CallSet}));
    }
    if ($hash{SomaticCallSet}) {
      @callers = (@callers, split(/\,/, $hash{SomaticCallSet}));
    }
  }
  $is_gs = 1;
  $is_gs = 0 if (grep(/hotspot/,@callers) && $maf > 0.1);
  $is_gs = 0 if ($maf < 0.05);
  $is_gs = 0 if ($maf < 0.1  && $vartype ne 'snp');
  $is_gs = 0 if ($hash{CallSetInconsistent} && $vartype ne 'snp');
  if ($id =~ m/COS/ && $cosmicsubj >= 5) {
    $is_gs = 0 if ($dp < 20 || $altct < 3);
  }elsif ($id =~ m/rs/){
      $is_gs = 0 if ($dp < 10 || $maf < 0.15);
  }else {
    $is_gs = 0 unless (scalar(@callers) > 1);
    $is_gs = 0 if ($dp < 20 || $altct < 8);
  }
  if ($is_gs) {
      print GS $line,"\n";
    }
  %chash = map {$_=>1} @callers;
  $vartype{$chrom}{$pos} = $vartype;
  if ($hash{PlatRef} =~ m/giab|platinum/) {
    if ($is_gs){
      $status{$chrom}{$pos}{'genomeseer'} = 'tp';
    }elsif($hash{PlatRef} =~ m/giab/ && $hash{PlatRef} =~ m/platinum/) {
      $status{$chrom}{$pos}{'genomeseer'} = 'fn' unless ($status{$chrom}{$pos}{'genomeseer'});
    }
    $status{$chrom}{$pos}{'union'} = 'tp';
    foreach (@allcallers) {
      if ($chash{$_}){
	$status{$chrom}{$pos}{$_} = 'tp';
      }elsif($hash{PlatRef} =~ m/giab/ && $hash{PlatRef} =~ m/platinum/) {
	$status{$chrom}{$pos}{$_} = 'fn' unless ($status{$chrom}{$pos}{'genomeseer'});
      }
    }
  }else {
    if ($is_gs){
      $status{$chrom}{$pos}{'genomeseer'} = 'fp';
      print FP $line,"\n";
    }
    $status{$chrom}{$pos}{'union'} = 'fp';
    foreach $caller (keys %chash) {
      $status{$chrom}{$pos}{$caller} = 'fp';
    }
  }
}
close FP;
close GS;
foreach $chrom (keys %status) {
  foreach $pos (keys %{$status{$chrom}}) {
    my $vartype = $vartype{$chrom}{$pos};
    foreach $call (keys %{$status{$chrom}{$pos}}) {
      if ($status{$chrom}{$pos}{$call} eq 'tp') {
	$tp{$vartype}{$call} ++;
      }elsif ($status{$chrom}{$pos}{$call} eq 'fp') {
	$fp{$vartype}{$call} ++;
      }elsif ($status{$chrom}{$pos}{$call} eq 'fn') {
	$fn{$vartype}{$call} ++;
      }
    }
  }
}

open OUT, ">$run\.snsp.txt" or die $!;
print OUT join("\t",'SampleID','Method','SNV TP','SNV FP','SNV FN','InDel TP','InDel FP','InDel FN',
	       'SNV SN','InDel SN','Combined SN','SNV PPV','InDel PPV','Combined PPV'),"\n";

foreach ((@allcallers,'genomeseer','union')) {
  $fp{'snp'}{$_} = 0 unless $fp{'snp'}{$_};
  $fp{'indel'}{$_} = 0 unless $fp{'indel'}{$_};
  $fn{'indel'}{$_} = 0 unless $fn{'indel'}{$_};
  next unless ($tp{'snp'}{$_});
  next unless ($tp{'indel'}{$_});
  $sn_snp = sprintf("%.1f",100*$tp{'snp'}{$_}/($tp{'snp'}{$_}+$fn{'snp'}{$_}));
  $sn_indel = sprintf("%.1f",100*$tp{'indel'}{$_}/($tp{'indel'}{$_}+$fn{'indel'}{$_}));
  $sn_total = sprintf("%.1f",100*($tp{'snp'}{$_}+$tp{'indel'}{$_})/($tp{'snp'}{$_}+$fn{'snp'}{$_}+$tp{'indel'}{$_}+$fn{'indel'}{$_}));
  $sp_total = sprintf("%.1f",100*($tp{'snp'}{$_}+$tp{'indel'}{$_})/($tp{'snp'}{$_}+$fp{'snp'}{$_}+$tp{'indel'}{$_}+$fp{'indel'}{$_}));
  $sp_snp =   sprintf("%.1f",100*$tp{'snp'}{$_}/($tp{'snp'}{$_}+$fp{'snp'}{$_}));
  $sp_indel = sprintf("%.1f",100*$tp{'indel'}{$_}/($tp{'indel'}{$_}+$fp{'indel'}{$_}));
  
  print OUT join("\t",$run,$_,$tp{'snp'}{$_},$fp{'snp'}{$_},$fn{'snp'}{$_},
		 $tp{'indel'}{$_},$fp{'indel'}{$_},$fn{'indel'}{$_},
		 $sn_snp,$sn_indel,$sn_total,$sp_snp,$sp_indel,$sp_total),"\n";
}
close OUT;
