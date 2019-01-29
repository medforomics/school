#!/usr/bin/perl -w
#integrate_datasets.pl

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %opt = ();
my $results = GetOptions (\%opt,'subject|s=s','rnaseqid|r=s','vcf|v=s','help|h');

open OUT, ">$opt{subject}\.itd.vcf" or die $!;

my @sampids;
open IN, "<$opt{vcf}" or die $!;
W1:while (my $line = <IN>) {
  chomp($line);
  if ($line =~ m/^#CHROM/) {
    my @header = split(/\t/,$line);
    ($chrom, $pos,$id,$ref,$alt,$score,
     $filter,$info,$format,@gtheader) = split(/\t/, $line);
    @sampids = @gtheader;
    push @sampids, $opt{rnaseqid};
    print OUT join("\t",$chrom,$pos,$id,$ref,$alt,$score,
		   $filter,$info,$format,@sampids),"\n";
    next;
  } elsif ($line =~ m/^#/) {
    print OUT $line,"\n";
    next;
  }
  my ($chrom, $pos,$id,$ref,$alt,$score,
      $filter,$annot,$format,@gts) = split(/\t/, $line);

  next if ($ref =~ m/\./ || $alt =~ m/\./ || $alt=~ m/,X/);
  push @gts, join(":",'.','.','.','.','.');
  print OUT join("\t",$chrom, $pos,$id,$ref,$alt,$score,$filter,$annot,
		 $format,@gts),"\n";
}
close IN;
close OUT;
