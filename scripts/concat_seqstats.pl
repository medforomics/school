#!/usr/bin/perl -w

open LIST, "<seqstatfiles.txt" or die $!;
while (my $sfile = <LIST>) {
    chomp($sfile);
    open IN, "$sfile" or die $!;
    while (my $line = <IN>) {
	chomp($line);
	my ($key,$value) = split(/\t/,$line);
	$columns{$key} = 1;
	$info{$sfile}{$key} = $value;
    }
}
open OUT, ">sequence.stats.txt" or die $!;
my @cols = sort {$a cmp $b} keys %columns;
print OUT join("\t",@cols),"\n";
foreach $dset (keys %info) {
    my @line;
    foreach $col (@cols) {
	push @line, $info{$dset}{$col};
    }
    print OUT join("\t",@line),"\n";
}
close OUT;
