#!/usr/bin/perl -w
#run_casava.pl

my $seqdatadir = shift @ARGV;

open IN, "<$seqdatadir\/RunInfo.xml.ori" or die $!;
open OUT, ">$seqdatadir\/RunInfo.xml" or die $!;
while (my $line = <IN>) {
    chomp($line);
    if ($line =~ m/(\s+)<Read Number="(\d+)"\s+NumCycles="(\d+)"/) {
	if ($2 eq 1) {
	    print OUT $line,"\n";
	} elsif ($2 eq 2) {
	    $numcycles = $3-8;
	    print OUT $1.qq{<Read Number="2" NumCycles="$numcycles" IsIndexedRead="Y" />},"\n";
	} elsif ($2 eq 3) {
	    $numcycles = $3+8;
	    print OUT $1.qq{<Read Number="3" NumCycles="$numcycles" IsIndexedRead="N" />},"\n";
	} 
    }else {
	print OUT $line,"\n";
    }
}

