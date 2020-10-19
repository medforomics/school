#!/usr/bin/perl -w
#validation_json2txt.pl

use JSON;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my %opt = ();
my $results = GetOptions (\%opt,'jsonfile|i=s');

my @vtypes = ('variants','cnvs','translocations');
$jsontxt = `cat $opt{jsonfile}`;
$jsonref  = decode_json($jsontxt);

foreach $case (@{$jsonref}) {
    open OUT, ">$case->{caseName}.reported.targets.txt" or die $!;
    foreach $tref (@{$case->{'variants'}}) {
	my %hash = %{$tref};
	print OUT join("\t",$hash{chrom},$hash{pos},$hash{reference},$hash{alt}),"\n";
    }
    close OUT;
}

