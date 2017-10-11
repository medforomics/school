#!/usr/bin/perl -w
#extract_greyzone.pl

my $vcffile = shift @ARGV;
my $prefix = $vcffile;
$prefix =~ s/\.vcf.gz//;
my $input = "$vcffile" or die $!;
open OUT, ">$prefix\.tumornormal.txt" or die $!;
print OUT join("\t",'CHROM','POS','ID','Gene','AminoAcid',
	       'TumorAF','TumorDepth','NormalAF','NormalDepth',
	       'RnaSeqAF','RnaSeqDepth'),"\n";

open IN, "gunzip -c $input|" or die $!;
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
  next unless ($filter eq 'PASS');
  my %hash = ();
  foreach $a (split(/;/,$annot)) {
    my ($key,$val) = split(/=/,$a);
    $hash{$key} = $val unless ($hash{$key});
  }
  $hash{NormalAF} = $hash{NormalMAF} if ($hash{NormalMAF});
  $hash{RnaSeqDP} = 0 unless $hash{RnaSeqDP};
  $hash{RnaSeqMAF} = 0 unless $hash{RnaSeqMAF};

  foreach $trx (split(/,/,$hash{ANN})) {
    my ($allele,$effect,$impact,$gene,$geneid,$feature,
	$featureid,$biotype,$rank,$codon,$aa,$pos_dna,$len_cdna,
	$cds_pos,$cds_len,$aapos,$aalen,$distance,$err) = split(/\|/,$trx);
    next unless ($impact =~ m/HIGH|MODERATE/);
    print OUT join("\t",$chrom,$pos,$id,$gene,$aa,$hash{AF},$hash{DP},$hash{NormalAF},$hash{NormalDP},$hash{RnaSeqMAF},$hash{RnaSeqDP}),"\n";
    next W1;
  }
}
