#!/usr/bin/perl -w
#integrate_datasets.pl

my $vcffile = shift @ARGV;
my $philipsvcf = $vcffile;
$philipsvcf =~ s/utsw.vcf.gz/philips.vcf/;
open IN, "gunzip -c $vcffile |" or die $!;
open OUT, ">$philipsvcf" or die $!;

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
  next if ($somline{$chrom}{$pos});
  my %hash = ();
  foreach $a (split(/;/,$annot)) {
    my ($key,$val) = split(/=/,$a);
    next if ($key =~ m/^Normal|RnaSeqAF|RnaSeqDP|RnaSeqAD|Germline/);
    $hash{$key} = $val unless ($hash{$key});
  }
  my @nannot;
  foreach $info (sort {$a cmp $b} keys %hash) {
    if (defined $hash{$info}) {
      push @nannot, $info."=".$hash{$info};
    }else {
      push @nannot, $info;
    }
  }
  $newannot = join(";",@nannot);
  print OUT join("\t",$chrom, $pos,$id,$ref,$alt,$score,
		 $filter,$newannot,$format,@gts),"\n";
}
