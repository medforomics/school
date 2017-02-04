#!/usr/bin/perl -w
#compile_gene_fusion_results.pl

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my $results = GetOptions (\%opt,'input|i=s','refdb|r=s','sname|n=s','help|h');

my $prefix = $opt{input};

open VCF, "<$prefix\.bk_sv.vcf" or die $!;
open BED, ">$prefix\.bed" or die $!;
my %bed;
while (my $line = <VCF>) {
  chomp($line);
  next if $line =~ m/#/;
  my ($chrom, $pos,$id,$ref,$alt,$score,
      $filter,$annot,$format,@gts) = split(/\t/, $line);
  $id =~ s/\//\|/g;
  my $evid = (split(/_/,$id))[1];
  my $end = (split(/_/,$id))[-1];
  print BED join("\t",'chr'.$chrom,$pos-1,$pos,$evid."_".$end),"\n";
}
close BED;
system(qq{bedtools intersect -header -wb -a $prefix\.bed -b $opt{refdb}\/gencode.exons.bed > exonoverlap_sv.txt});
system(qq{bedtools intersect -header -v -wb -a $prefix\.bed -b $opt{refdb}\/gencode.exons.bed | bedtools intersect -header -wb -a stdin -b $opt{refdb}\/gencode.genes.chr.bed > geneoverlap_sv.txt});


open IN, "<$prefix\.summary.tsv" or die $!;
my $head = <IN>;
chomp($head);
my @header = split(/\t/,$head);
my %readct;
while (my $line = <IN>) {
  chomp($line);
  my @row = split(/\t/,$line);
  my %hash;
  foreach my $j (0..$#row) {
    $hash{$header[$j]} = $row[$j];
  }
  $readct{$hash{Fusion_Candidate}} = $hash{EN_RNA}+$hash{SP_RNA}+$hash{EN_DNA_T}+$hash{SP_DNA_T};
}

open IN, "<exonoverlap_sv.txt" or die $!;
while (my $line = <IN>) {
  chomp($line);
  next if ($line =~ m/#/);
  my ($chrom,$start,$end,$id,$chr,$start2,$end2,$gene) = split(/\t/, $line);
  my $evid = (split(/_/,$id))[0];
  my %done;
  foreach $exon (split(/,/,$gene)) {
    my ($symbol,$trxid,$exonnum) = split(/\|/,$exon);
    next if ($done{$symbol}{$exonnum});
    $done{$symbol}{$exonnum} = 1;
    push @{$fusions{$evid}{$symbol}}, join("\t",join("_",'gf',$evid,$opt{sname}),
					   $chrom,$start,$end,
					   'gene_fusion',$readct{$evid},
					   $symbol,$exonnum,$opt{sname});
  }
}
open IN, "<geneoverlap_sv.txt" or die $!;
while (my $line = <IN>) {
  chomp($line);
  next if ($line =~ m/#/);
  my ($chrom,$start,$end,$id,$chr,$start2,$end2,$gene) = split(/\t/, $line);
  my $evid = (split(/_/,$id))[0];
  my ($symbol,$ensid) = split(/\|/,$gene);
  push @{$fusions{$evid}{$symbol}}, join("\t",join("_",'gf',$evid,$opt{sname}),
					 $chrom,$start,$end,
					 'gene_fusion',$readct{$evid},
					 $symbol,'',$opt{sname});
}
open OUT, ">$prefix\.genefusions.txt" or die $!;
foreach $evid (sort {$a <=> $b} keys %fusions) {
  my @genes = keys %{$fusions{$evid}};
  if (scalar(@genes) > 1) {
    foreach $sym (@genes) {
      print OUT join("\n", @{$fusions{$evid}{$sym}}),"\n";
    }
  }
}

