#!/usr/bin/perl -w
#integrate_datasets.pl

my ($somaticvcf) = @ARGV;
my ($tumorid,$normalid);

my $outfile = $somaticvcf;
$outfile =~ s/annot.vcf.gz/somatic.vcf/;

die if ($somaticvcf eq $outfile);

open OUT, ">$outfile" or die $!;
open IN, "gunzip -c $somaticvcf |" or die $!;
W1:while (my $line = <IN>) {
  chomp($line);
  if ($line =~ m/^#CHROM/) {
    my @header = split(/\t/,$line);
    ($chrom, $pos,$id,$ref,$alt,$score,
     $filter,$info,$format,@gtheader) = split(/\t/, $line);
    print OUT join("\t",$chrom,$pos,$id,$ref,$alt,$score,
		   $filter,$info,$format,@gtheader),"\n";

    foreach $id (@gtheader) {
      if ($id =~ m/_T_/) {
	$tumorid = $id;
      }else {
	$normalid = $id;
      }
    }
    next;
  } elsif ($line =~ m/^#/) {
    $line =~ s/CallSet/SomaticCallSet/;
    print OUT $line,"\n";
    next;
  }
  my ($chrom, $pos,$id,$ref,$alt,$score,
      $filter,$annot,$format,@gts) = split(/\t/, $line);
  next if ($ref =~ m/\./ || $alt =~ m/\./ || $alt=~ m/,X/);
  $annot =~ s/CallSet/SomaticCallSet/;
  my %hash = ();
  foreach $a (split(/;/,$annot)) {
    my ($key,$val) = split(/=/,$a);
    $hash{$key} = $val unless ($hash{$key});
  }
  my %fail;
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
      $hash{$deschead[$k]} = $ainfo[$k] if ($subjid eq $tumorid);
    }
    $gtinfo{$subjid}{DP} = (split(/,/,$gtinfo{$subjid}{DP}))[0] if ($gtinfo{$subjid}{DP});
    next F1 unless ($gtinfo{$subjid}{DP} && $gtinfo{$subjid}{DP} ne '.' && $gtinfo{$subjid}{DP} >= 5);
    my @altct = split(/,/,$gtinfo{$subjid}{AO});
    foreach  my $act (@altct) {
      next if ($act eq '.');
      push @mutallfreq, sprintf("%.4f",$act/$gtinfo{$subjid}{DP});
    }
    $gtinfo{$subjid}{MAF} = \@mutallfreq;
  }
  next unless ($gtinfo{$tumorid}{DP} && $gtinfo{$tumorid}{DP} >= 20);
  @tumormaf = @{$gtinfo{$tumorid}{MAF}};
  @tumoraltct = split(/,/,$gtinfo{$tumorid}{AO});
  if (exists $hash{INDEL}) {
      $hash{TYPE} = 'indel';
  }
  $hash{TYPE} = 'ambi' unless ($hash{"TYPE"});
  next if ($tumoraltct[0] eq '.');
  $hash{AF} = join(",",@tumormaf);
  my @callers = split(/,/,$hash{SomaticCallSet});
  $hash{SS} = 5;
  delete $hash{SOMATIC};
  next unless ($gtinfo{$normalid} && exists $gtinfo{$normalid}{MAF});
  @normalmaf = @{$gtinfo{$normalid}{MAF}};
  $hash{NormalAF} = $normalmaf[0];
  $hash{NormalDP} = $gtinfo{$normalid}{DP};
  if ($normalmaf[0] >= 0.30) {
    next;
  }elsif ($tumormaf[0] < 0.05 && ($normalmaf[0] > 0.01 || $normalmaf[0]*5 > $tumormaf[0])) {
    next;
  }elsif ($normalmaf[0] >= 0.05 || $normalmaf[0]*5 > $tumormaf[0]) {
    next;
  }else {
    $hash{SS} = 2;
  }
  my $newformat = 'GT:DP:AD:AO:RO';
  my @newgt;
  foreach $sample (@gtheader) {
    my @gtdata;
      foreach $gt (split(/:/,$newformat)) {
	  $gtinfo{$sample}{$gt} = '.' unless (exists $gtinfo{$sample}{$gt});
	  push @gtdata, $gtinfo{$sample}{$gt};
      }
      push @newgt, join(":",@gtdata);
  }
  print OUT join("\t",$chrom, $pos,$id,$ref,$alt,$score,$filter,$annot,
		 $newformat,@newgt),"\n";
}
close IN;
