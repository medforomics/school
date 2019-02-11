#!/usr/bin/perl -w
#use strict;

open OM, "</project/shared/bicf_workflow_ref/human/GRCh38/panel1385.genelist.txt" or die $!;
while (my $line = <OM>) {
  chomp($line);
  $keep{$line} = 1;
}

my $vcffile = shift @ARGV;
my $prefix = $vcffile;
$prefix =~ s/\.vcf.gz//;
my $input = "$vcffile" or die $!;
open OUT, ">$prefix\.tumor.vcf" or die $!;
open IN, "gunzip -c $input|" or die $!;
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
  #unless ($exacaf eq '' || $exacaf <= 0.01) {
  #  $fail{'COMMON'} = 1;
  #}
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
  foreach  my $i (0..$#acts) {
    if ($acts[$i] > 2) {
      push @newalts, $altnts[$i];
      push @altct, $acts[$i];
      $totalaltct += $acts[$i];
    }
  }
  next unless ($altct[0] && $altct[0] > 2);
  if ($hash{DP} =~ m/,/) {
    $hash{DP} = $totalaltct+$hash{RO};
  }
  next unless ($hash{DP});
  my @mutallfreq;
  foreach  my $act (@altct) {
    push @mutallfreq, sprintf("%.4f",$act/$hash{DP});
  }
  my @sortao = sort {$b <=> $a} @altct;
  $hash{AF} = join(",",@mutallfreq);
  if ($hash{DP} < 20) {
    $fail{'LowDepth'} = 1;
  }
  my @callers = split(/,/,$hash{CallSet});
  if ($somval{$chrom}{$pos}) {
    push @callers, split(/,/,$somval{$chrom}{$pos});
  }
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
    next unless ($impact =~ m/HIGH|MODERATE/);
    next unless $keep{$gene};
    push @aa, $aa if ($aa ne '');
    $keepforvcf = $gene;
  }
  next unless $keepforvcf;
  my @fail = keys %fail;
  if (scalar(@fail) < 1) {
    $filter = 'PASS';
    print OUT join("\t",$chrom,$pos,$id,$ref,join(",",@newalts),$score,
		   $filter,$annot,$format,$allele_info),"\n";
  }else {
    next;
  }
}
