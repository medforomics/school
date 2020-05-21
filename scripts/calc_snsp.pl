#!/usr/bin/perl -w
#calc_snp.pl

my %tp;
my %fp;
my %fn;

my $input = shift @ARGV;
my $fn = 'fn.vcf';
my $run = $input;
$run =~ s/\.eval.vcf//;

my @allcallers = ('strelka2','fb','platypus','gatk','genomeseer','union');
foreach ((@allcallers)) {
  $fn{'snp'}{$_} = 0;
  $fn{'indel'}{$_} = 0;
}
open FN, "<$fn" or die $!;
while (my $line = <FN>) {
  chomp($line);
  my ($chr,$pos,$id,$ref,$alt,@others) = split(/\t/,$line);
  $vartype = 'snp';
  if (length($ref) > 1 || length($alt) > 1) {
    $vartype = 'indel';
  }
  next if ($chr eq 'chr6' && $pos >= 28510120 && $pos <= 33480577);
  foreach (@allcallers) {
    $fn{$vartype}{$_} ++;
  }
}

my %done;
my %reason;
open IN, "<$input" or die $!;
while (my $line = <IN>) {
  chomp($line);
  if ($line =~ m/^#/) {
    next;
  }
  my ($chrom, $pos,$id,$ref,$alt,$score,
      $filter,$annot,$format,@gts) = split(/\t/, $line);
  next if ($ref =~ m/\./ || $alt =~ m/\./ || $alt=~ m/,X/);
  next if ($chrom eq 'chr6' && $pos >= 28510120 && $pos <= 33480577);
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
  if ($maf < 0.05 || ($maf < 0.10  && $vartype ne 'snp')) {
      $reason{$chrom}{$pos}{$ref}{$alt} = 'LowAF';
      $is_gs = 0;
  }elsif ($hash{CallSetInconsistent} && $vartype ne 'snp') {
      $reason{$chrom}{$pos}{$ref}{$alt} = 'Inconsistent';
      $is_gs = 0;
  }elsif ($hash{RepeatType} && $hash{RepeatType} =~ m/Simple_repeat/ && $maf < 0.15) {
      $reason{$chrom}{$pos}{$ref}{$alt} = 'InRepeat';
      $is_gs = 0;
  } elsif (($hash{FS} && $hash{FS} > 60) || $filter =~ m/strandBias/i || $hash{strandBias} || (($hash{SAP} && $hash{SAP} > 20) && ((exists $hash{SAF} && $hash{SAF}< 1 & $hash{SRF} > 10) || (exists $hash{SAR} && $hash{SAR}< 1 & $hash{SRR} > 10))))
  {
      $reason{$chrom}{$pos}{$ref}{$alt} = 'StrandBias';
      $is_gs = 0;
  } elsif ($id =~ m/COS/ && $cosmicsubj >= 5) {
      if ($dp < 20 || $altct < 3) {
	  $reason{$chrom}{$pos}{$ref}{$alt} = 'LowCoverage';
	  $is_gs = 0;
      }
  } elsif  (scalar(@callers) < 2) {
      $reason{$chrom}{$pos}{$ref}{$alt} = 'OneCaller';
      $is_gs = 0;
  } elsif ($dp < 20 || $altct < 8) {
      $reason{$chrom}{$pos}{$ref}{$alt} = 'LowCoverage';
      $is_gs = 0;
  }
  my %sinfo;
  %chash = map {$_=>1} @callers;
  if ($is_gs) {
    $chash{'genomeseer'} = 1;
  }
  $chash{'union'} = 1;
  $vartype{$chrom}{$pos}{$ref}{$alt} = $vartype;
  $refpos = 0;
  if ($hash{platforms}) {
      $refpos ++;
  }if ($hash{MTD}) {
      $refpos ++;
  }
  my $chrompos= join(":",$chrom,$pos,$ref,$alt);
  my $callers = join(":",keys %chash);
  $status{$chrompos}{$callers} = $refpos;
}
my %cat;
foreach $loci (keys %status) {
  my @csets = keys %{$status{$loci}};
  my ($chrom,$pos,$ref,$alt) = split(/:/,$loci);
  my $trueset = 0;
  my $genomeseer = 0;
  my %callers;
  foreach $cset (@csets) {
    $repos = $status{$loci}{$cset};
    @callers = split(/:/,$cset);
    if ($trueset) {	
      $trueset = $repos if ($repos > $trueset); 
    }else {
      $trueset = $repos;
    }
    foreach (@callers) {
      $callers{$_} = 1;
    }
  }
  if ($trueset > 0) {
    foreach (@allcallers) {
      if ($callers{$_}){
	$cat{$chrom}{$pos}{$ref}{$alt}{$_} = 'tp';
      } else {
	$cat{$chrom}{$pos}{$ref}{$alt}{$_} = 'fn' unless ($trueset < 2);
	print join(":","FN",$chrom,$pos,$ref,$alt,$reason{$chrom}{$pos}{$ref}{$alt}),"\n" if ($_ eq 'genomeseer');
      }
    }
  }else {
    foreach (@allcallers) {
      if ($callers{$_}){
	$cat{$chrom}{$pos}{$ref}{$alt}{$_} = 'fp';
	print join(":","FP",$chrom,$pos,$ref,$alt,$vartype{$chrom}{$pos}{$ref}{$alt}),"\n" if ($_ eq 'genomeseer');
      }
    }
  }
}
foreach $chrom (keys %cat) {
  foreach $pos (keys %{$cat{$chrom}}) {
	foreach $ref (keys %{$cat{$chrom}{$pos}}){
		foreach $alt (keys %{$cat{$chrom}{$pos}{$ref}}){
		    my $vartype = $vartype{$chrom}{$pos}{$ref}{$alt};
		    foreach $call (keys %{$cat{$chrom}{$pos}{$ref}{$alt}}) {
  			    if ($cat{$chrom}{$pos}{$ref}{$alt}{$call} eq 'tp') {
					$tp{$vartype}{$call} ++;
	      		}elsif ($cat{$chrom}{$pos}{$ref}{$alt}{$call} eq 'fp') {
					$fp{$vartype}{$call} ++;
				}elsif ($cat{$chrom}{$pos}{$ref}{$alt}{$call} eq 'fn') {
					$fn{$vartype}{$call} ++;
				}
			}
    	}
    }
  }
}

open OUT, ">$run\.snsp.txt" or die $!;
print OUT join("\t",'SampleID','Method','SNV TP','SNV FP','SNV FN','InDel TP','InDel FP','InDel FN',
	       'SNV SN','InDel SN','Combined SN','SNV PPV','InDel PPV','Combined PPV'),"\n";

foreach (@allcallers) {
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
