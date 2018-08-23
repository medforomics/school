#!/usr/bin/perl -w
#integrate_datasets.pl

#module load vcftools/0.1.14 samtools/1.6 bedtools/2.26.0 
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %opt = ();
my $results = GetOptions (\%opt,'subject|s=s','id|i=s','vcf|v=s','help|h');

open OUT, ">$opt{subject}\.all.vcf" or die $!;
open PASS, ">$opt{subject}\.pass.vcf" or die $!;

my @sampids;
open IN, "gunzip -c $opt{vcf} |" or die $!;
W1:while (my $line = <IN>) {
  chomp($line);
  if ($line =~ m/^#CHROM/) {
    print OUT qq{##INFO=<ID=RnaSeqAF,Number=A,Type=Float,Description="RNASeq Allele Frequency">\n};
    print OUT qq{##INFO=<ID=RnaSeqDP,Number=1,Type=Integer,Description="RNASeq read depth">\n};
    print PASS qq{##INFO=<ID=RnaSeqAF,Number=A,Type=Float,Description="RNASeq Allele Frequency">\n};
    print PASS qq{##INFO=<ID=RnaSeqDP,Number=1,Type=Integer,Description="RNASeq read depth">\n};
    my @header = split(/\t/,$line);
    ($chrom, $pos,$id,$ref,$alt,$score,
     $filter,$info,$format,@gtheader) = split(/\t/, $line);
    unless ($opt{id}) {
	$opt{id} = $gtheader[0];
    }
    @sampids = ($opt{id});
    push @sampids, $rnaseqid if ($rnaseqid);
    print OUT join("\t",$chrom,$pos,$id,$ref,$alt,$score,
		   $filter,$info,$format,@sampids),"\n";
    print PASS join("\t",$chrom,$pos,$id,$ref,$alt,$score,
		    $filter,$info,$format,@sampids),"\n";
    next;
  } elsif ($line =~ m/^#/) {
    print OUT $line,"\n";
    print PASS $line,"\n";
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
  $hash{'HG38Loci'} = join(":",$chrom,$pos);
  my %fail;
  $fail{'UTSWBlacklist'} = 1 if ($hash{UTSWBlacklist});
  my @exacaf;
  my $exacaf;
  if ($hash{AF_POPMAX}) {
    foreach (split(/,/,$hash{AF_POPMAX})) {
      push @exacaf, $_ if ($_ ne '.');
    }
    @exacaf = sort {$b <=> $a} @exacaf;
    $exacaf = $exacaf[0] if ($exacaf[0]);
  }elsif ($hash{dbNSFP_ExAC_Adj_AF}) {
    foreach (split(/,/,$hash{dbNSFP_ExAC_Adj_AF})) {
      push @exacaf, $_ if ($_ ne '.');
    }
    @exacaf = sort {$b <=> $a} @exacaf;
    $exacaf = $exacaf[0] if ($exacaf[0]);
  }elsif ($hash{AC_POPMAX} && $hash{AN_POPMAX}) {
    my @exacs = split(/,/,$hash{AC_POPMAX});
    my $ac = 0;
    foreach $val (@exacs) {
      $ac += $val if ($val =~ m/^\d+$/);
    }
    my @exans = split(/,/,$hash{AN_POPMAX});
    my $an = 0;
    foreach $val (@exans) {
      $an += $val if ($val =~ m/^\d+$/);
    }
    $exacaf = sprintf("%.4f",$ac/$an) if ($ac > 0 && $an > 10);
  }
  $fail{'COMMON'} = 1 if ($exacaf && $exacaf > 0.01);
  next if ($exacaf && $exacaf > 0.2);
  $fail{'StrandBias'} = 1 if (($hash{FS} && $hash{FS} > 60) || $filter =~ m/strandBias/i);
  my $cosmicsubj = 0;
  if ($hash{CNT}) {
    my @cosmicct = split(/,/,$hash{CNT}); 
    foreach $val (@cosmicct) {
      $cosmicsubj += $val if ($val =~ m/^\d+$/);
    }
  }
  my %gtinfo = ();
  my @deschead = split(/:/,$format);
 F1:foreach my $k (0..$#gtheader) {
    my $subjid = $gtheader[$k];
    my $allele_info = $gts[$k];
    my @ainfo = split(/:/, $allele_info);
    my @mutallfreq = ();
    foreach my $k (0..$#ainfo) {
      $gtinfo{$subjid}{$deschead[$k]} = $ainfo[$k];
      $hash{$deschead[$k]} = $ainfo[$k] if ($subjid eq $opt{id});
    }
    $gtinfo{$subjid}{DP} = (split(/,/,$gtinfo{$subjid}{DP}))[0] if ($gtinfo{$subjid}{DP});
    next F1 unless ($gtinfo{$subjid}{DP} && $gtinfo{$subjid}{DP} ne '.' && $gtinfo{$subjid}{DP} >= 1);
    my @altct = split(/,/,$gtinfo{$subjid}{AD});
    my $refct = shift @altct;
    @altct2 = split(/,/,$gtinfo{$subjid}{AO});
    if (scalar(@altct) ne scalar(@altct2)) {
	warn "Inconsistent Allele counts @ $chrom,$pos,$alt,$gtinfo{$subjid}{AD},$gtinfo{$subjid}{AO}\n";
    }
    my $total = $refct;
    foreach  my $act (@altct) {
	next if ($act eq '.');
	$total += $act;
	push @mutallfreq, sprintf("%.4f",$act/$gtinfo{$subjid}{DP});
    }
    $gtinfo{$subjid}{MAF} = \@mutallfreq;
  }
  next unless ($gtinfo{$opt{id}}{DP} && $gtinfo{$opt{id}}{DP} ne '.' && $gtinfo{$opt{id}}{DP} >= 20);
  unless ($gtinfo{$opt{id}}{AO} =~ m/\d+/ && $gtinfo{$opt{id}}{AD} =~ m/,/) {
      warn "Missing Alt:$line\n";
  }
  @tumormaf = @{$gtinfo{$opt{id}}{MAF}};
  @tumoraltct = split(/,/,$gtinfo{$opt{id}}{AO});
  
  if (exists $hash{INDEL}) {
      $hash{TYPE} = 'indel';
  }
  $hash{TYPE} = 'ambi' unless ($hash{"TYPE"});
  next if ($tumoraltct[0] eq '.');
  $hash{AF} = join(",",@tumormaf);
  my @callers;
  if ($hash{CallSet} && $hash{CallSet} =~ m/\// || $hash{SomaticCallSet} && $hash{SomaticCallSet} =~ m/\//) {
      my @callinfo ;
      @callinfo = split(/\|/, $hash{CallSet}) if ($hash{CallSet});
      if ($hash{SomaticCallSet}) {
	  @callinfo = (@callinfo, split(/\|/, $hash{SomaticCallSet}));
      }
      foreach $cinfo (@callinfo) {
	  my ($caller, $alt, @samafinfo) = split(/\//,$cinfo);
	  push @callers, $caller;
      }
      $hash{CallSet} = join(",",@callinfo);
      $hash{CallSet} =~ s/\//\|/g;
  } elsif ($hash{CallSet} && $hash{CallSet} =~ m/\|/ || $hash{SomaticCallSet} && $hash{SomaticCallSet} =~ m/\|/) {
      my @callinfo ;
      @callinfo = split(/,/, $hash{CallSet}) if ($hash{CallSet});
      if ($hash{SomaticCallSet}) {
	  @callinfo = (@callinfo, split(/,/, $hash{SomaticCallSet}));
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
      $hash{CallSet} = join(",",@callers);
  }
  if ($id =~ m/COS/ && $cosmicsubj >= 5) {
      $fail{'LowAltCt'} = 1 if ($tumoraltct[0] < 3);
      $fail{'LowMAF'} = 1 if ($tumormaf[0] < 0.05);
      $fail{'LowMAF'} = 1 if ($tumormaf[0] < 0.1 && $hash{TYPE} ne 'snp');
  }else {
    $fail{'OneCaller'} = 1 if (scalar(@callers) < 2);
    $fail{'LowAltCt'} = 1 if ($tumoraltct[0] < 8);
    $fail{'LowMAF'} = 1 if ($tumormaf[0] < 0.05);
    $fail{'LowMAF'} = 1 if ($tumormaf[0] < 0.10 && $hash{TYPE} ne 'snp');
  }
  delete $hash{SOMATIC};
  $hash{SS} = 5  unless ($hash{SS});
  if ($rnaval{$chrom}{$pos}) {
    $gtinfo{$rnaseqid} ={GT=>'.',DP=>'.',AO=>'.',AD=>'.',RO=>'.'};
    my ($rnahashref,$rnadp) = @{$rnaval{$chrom}{$pos}};
    if ($rnadp > 10) {
      my %rnantct = %{$rnahashref};
      my @altcts;
      my $totalaltct =0;
      foreach $altnt (split(/,/,$alt)) {
	my $ct = $rnantct{$altnt};
	$ct = 0 unless ($ct);
	$totalaltct += $ct;
	      push @altcts, $ct;
      }
      $hash{RnaSeqDP} = $rnadp;
      $hash{RnaSeqAF} = sprintf("%.4f",$altcts[0]/$rnadp);
      $gtinfo{$rnaseqid}{RO} = $rnadp - $totalaltct;
      $gtinfo{$rnaseqid}{AO} = join(",",@altcts);
      $gtinfo{$rnaseqid}{GT} = '.';
      $gtinfo{$rnaseqid}{DP} = $rnadp;
      $gtinfo{$rnaseqid}{AD} = join(",",$gtinfo{$rnaseqid}{RO},@altcts);
    }
  }
  if ($rnaseqid) {
      if ($gtinfo{$rnaseqid}{AO} && $gtinfo{$rnaseqid}{AO} =~ m/^\d+/ &&
	(split(/,/,$gtinfo{$rnaseqid}{AO}))[0] > 2) {
      $hash{RnaSeqValidation} = 1
    }
  }
  my $newformat = 'GT:DP:AD:AO:RO';
  my @newgt;
  foreach $sample (@sampids) {
    my @gtdata;
    foreach $gt (split(/:/,$newformat)) {
      $gtinfo{$sample}{$gt} = '.' unless (exists $gtinfo{$sample}{$gt});
      push @gtdata, $gtinfo{$sample}{$gt};
    }
    push @newgt, join(":",@gtdata);
  }
  next unless ($hash{ANN});
  my $cancergene = 0;
  my $keepforvcf;
  foreach $trx (split(/,/,$hash{ANN})) {
    my ($allele,$effect,$impact,$gene,$geneid,$feature,
	$featureid,$biotype,$rank,$codon,$aa,$pos_dna,$len_cdna,
	$cds_pos,$cds_len,$aapos,$aalen,$distance,$err) = split(/\|/,$trx);
    next unless ($impact =~ m/HIGH|MODERATE/ || $effect =~ /splice/i);
    next if($effect eq 'sequence_feature');
    $keepforvcf = $gene;
    #$cancergene = 1 if ($cgenelist{$gene});
  }
  next unless $keepforvcf;
  my @fail = sort {$a cmp $b} keys %fail;
  if (scalar(@fail) == 0) {
    $filter = 'PASS';
  }else {
    $filter = join(";", 'FailedQC',@fail);
  }
  my @nannot;
  foreach $info (sort {$a cmp $b} keys %hash) {
    if (defined $hash{$info}) {
      push @nannot, $info."=".$hash{$info};
    }else {
      push @nannot, $info;
    }
  }
  my $newannot = join(";",@nannot);
  print PASS join("\t",$chrom, $pos,$id,$ref,$alt,$score,$filter,$newannot,
		  $newformat,@newgt),"\n" if ($filter eq 'PASS');
  print OUT join("\t",$chrom, $pos,$id,$ref,$alt,$score,$filter,$newannot,
		 $newformat,@newgt),"\n" if ($filter eq 'PASS' || $id =~ m/COS/ || $cancergene || $filter eq 'FailedQC;COMMON');
}
close IN;
