#!/usr/bin/perl -w
#integrate_datasets.pl

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %opt = ();
my $results = GetOptions (\%opt,'outfile|o=s','rnaid|r=s',
			  'normal|n=s','vcf|v=s','help|h');

open OUT, ">$opt{outfile}" or die $!;

my @sampids;
my @addedgt;
open IN, "gunzip -c $opt{vcf} |" or die $!;
W1:while (my $line = <IN>) {
  chomp($line);
  if ($line =~ m/^#CHROM/) {
    my @header = split(/\t/,$line);
    ($chrom, $pos,$id,$ref,$alt,$score,
     $filter,$info,$format,@gtheader) = split(/\t/, $line);
    %sampids = map {$_=>1} @gtheader;
    unless ($sampids{$opt{normal}}) {
      push @gtheader, $opt{normal};
      push @addedgt, join(":",'.','.','.','.','.');
    }
    unless ($sampids{$opt{rnaid}}) {
      push @gtheader, $opt{rnaid};
      push @addedgt, join(":",'.','.','.','.','.');
    }
    print OUT join("\t",$chrom,$pos,$id,$ref,$alt,$score,
		   $filter,$info,$format,@gtheader),"\n";
    next;
  } elsif ($line =~ m/^#/) {
    print OUT $line,"\n";
    next;
  }
  my ($chrom, $pos,$id,$ref,$alt,$score,
      $filter,$annot,$format,@gts) = split(/\t/, $line);
  next if ($ref =~ m/\./ || $alt =~ m/\./ || $alt=~ m/,X/);
  push @gts, @addedgt;
  print OUT join("\t",$chrom, $pos,$id,$ref,$alt,$score,$filter,$annot,
		 $format,@gts),"\n";
}
close IN;
close OUT;
