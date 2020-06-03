#!/usr/bin/perl -w
#make_delly_sample.pl

my ($tumor,$normal) = @ARGV;
open OUT, ">samples.tsv" or die $!;
print OUT join("\t",$tumor,'tumor'),"\n";
print OUT join("\t",$normal,'control'),"\n";
close OUT;
