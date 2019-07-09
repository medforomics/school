#!/usr/bin/perl 
#parse_pindel

my $tumorid = shift @ARGV;
my $vcf = shift @ARGV;

open DUP, ">dupcnv.txt" or die $!;

open VCF, "gunzip -c $vcf|" or die $!;
while (my $line = <VCF>) {
  chomp($line);
  if ($line =~ m/#/) {
    if ($line =~ m/#CHROM/) {
      my @header = split(/\t/,$line);
      ($chrom, $pos,$id,$ref,$alt,$score,
       $filter,$info,$format,@subjacc) = split(/\t/, $line);
    }
    next;
  }
  my ($chrom, $pos,$id,$ref,$alt,$score,
      $filter,$annot,$format,@gts) = split(/\t/, $line);
  my %hash = ();
  foreach $a (split(/;/,$annot)) {
    my ($key,$val) = split(/=/,$a);
    $hash{$key} = $val;
  }
  my @deschead = split(/:/,$format);
  $hash{'END'} = $pos+1 unless $hash{'END'};
  my ($allele,$effect,$impact,$gene,$geneid,$feature,
      $featureid,$biotype,$rank,$codon,$aa,$pos_dna,$len_cdna,
      $cds_pos,$cds_len,$aapos,$aalen,$distance,$err);
  next unless $hash{ANN};
 F1:foreach $trx (split(/,/,$hash{ANN})) {
     ($allele,$effect,$impact,$gene,$geneid,$feature,
      $featureid,$biotype,$rank,$codon,$aa,$pos_dna,$len_cdna,
      $cds_pos,$cds_len,$aapos,$aalen,$distance,$err) = split(/\|/,$trx);
     next unless ($impact && $impact =~ m/HIGH|MODERATE/);
     next if ($effect eq 'sequence_feature');
     last F1;
 }
  next unless ($impact && $impact =~ m/HIGH|MODERATE/);
  next unless ($gene);
  next unless ($rank);
  next if ($gene eq 'FCGBP');
  print DUP join("\t",$gene,$chrom,$pos,$hash{END},'ITD',3,$hash{AF},$rank),"\n";
}
close VCF;
