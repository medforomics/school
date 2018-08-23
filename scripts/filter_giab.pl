#!/usr/bin/perl -w
#integrate_datasets.pl

my ($tumorid) = @ARGV;

open OUT, ">$tumorid\.PASS.vcf" or die $!;
open IN, "gunzip -c $tumorid\.annot.vcf.gz|" or die $!;
my %done;
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
  if ($hash{REFREP} && $hash{REFREP} =~ m/^\d$/) {
      $numreps= (split(/,/,$hash{REFREP}))[0];
      next if $numreps > 1;
  }
  my $vartype = 'snp';
  if (length($alt) == length($ref) && length($ref) > 1) {
      $vartype = 'indel';
  }
  my %fail;
  $fail{'StrandBias'} = 1 if (($hash{FS} && $hash{FS} > 60) || $filter =~ m/strandBias/i);
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
    if ($acts[$i] ne '.' && $acts[$i] > 2) {
      push @newalts, $altnts[$i];
      push @altct, $acts[$i];
      $totalaltct += $acts[$i];
    }
  }
  next unless ($altct[0] && $hash{DP} );
  next if ($hash{DP} < 20 || $altct[0] < 3);
  if ($hash{DP} =~ m/,/) {
      $hash{DP} = $totalaltct+$hash{RO};
  }
  $hash{AD} = join(",",$hash{RO},@altct);
  my @mutallfreq;
  foreach  my $act (@altct) {
      push @mutallfreq, sprintf("%.4f",$act/$hash{DP});
  }
  my @sortao = sort {$b <=> $a} @altct;
  $hash{AF} = join(",",@mutallfreq);
  next if ($vartype eq 'snp' && $mutallfreq[0] < 0.05);
  if ($vartype ne 'snp') {
      next if ( $mutallfreq[0] < 0.1)
  }
  my @callers;
  if ($hash{CallSet} && $hash{CallSet} =~ m/\|/ || $hash{SomaticCallSet} && $hash{SomaticCallSet} =~ m/\|/) {
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
  }else {
      if ($hash{CallSet}) {
	  @callers = (@callers, split(/\,/, $hash{CallSet}));
      }
      if ($hash{SomaticCallSet}) {
	  @callers = (@callers, split(/\,/, $hash{SomaticCallSet}));
      }
  }
  if ((grep(/hotspot/,@callers) || $id =~ m/COS/) && $cosmicsubj >= 5) {
      next if ($altct[0] < 3);
  }else {
      next if ($altct[0] < 8);
  }
  my @aa;
  next unless ($hash{ANN});
  foreach $trx (split(/,/,$hash{ANN})) {
      my ($allele,$effect,$impact,$gene,$geneid,$feature,
	  $featureid,$biotype,$rank,$codon,$aa,$pos_dna,$len_cdna,
	  $cds_pos,$cds_len,$aapos,$aalen,$distance,$err) = split(/\|/,$trx);
      push @aa, $aa if ($aa ne '');
  }

  next if (scalar(@fail) > 1);
  $filter = 'PASS';

  my @nannot;
  foreach $info (sort {$a cmp $b} keys %hash) {
      if (defined $hash{$info}) {
	  push @nannot, $info."=".$hash{$info};
      }else {
	  push @nannot, $info;
      }
  }
  $newannot = join(";",@nannot);
  print OUT join("\t",$chrom, $pos,$id,$ref,join(",",@newalts),$score,
		 $filter,$newannot,$format,$allele_info),"\n";
}
