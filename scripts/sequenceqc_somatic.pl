#!/usr/bin/perl -w
#uploadqc.pl

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %opt = ();
my $results = GetOptions (\%opt,'input|i=s','output|o=s','help|h');

open MATE, "<$opt{input}" or die $!;

while (my $line = <MATE>) {
    chomp($line);
    my ($sam1,$pf,$sam2,$corr,$depth) = split(/\t/,$line);
    $sam1 = (split(/\./,$sam1))[0];
    $sam2 = (split(/\./,$sam2))[0];
    open OUT, ">$opt{output}" or die $!;
    my $status= 'PASS';	
    $status='FAIL' if($pf eq 'unmatched');
    print OUT join("\n","Sample_1\t".$sam1,"Sample_2\t".$sam2,"Correlation\t".$corr,
		   "Depth\t".$depth,"Status\t".$status),"\n";
}
