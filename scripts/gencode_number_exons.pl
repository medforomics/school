#!/usr/bin/perl -w
#parse_gencode.pl


my $keep = shift @ARGV;
my %inc;
if ($keep) {
    open KEEP, "<$keep" or die $!;
    while (my $line = <KEEP>) {
	chomp($line);
	$inc{$line} = 1;
    }
}

open GCODE, "<gencode.exon.bed" or die $!;
open OUT, ">gencode.exon_numbered.bed" or die $!;
my %hash;
my %strand;
while (my $line = <GCODE>) {
    chomp($line);
    next if ($line =~ m/^#/);
    my ($chrom,$start,$end,$info) = split(/\t/,$line);
    my %keep_region;
    foreach $annot (split(/,/,$info)) {
	my ($gene,$transname,$enum,$strand) = split(/\|/,$annot);
	$strand{$gene} = $strand;
	if ($keep) {
	    next unless $inc{$gene};
	}
	$keep_region{$gene} = 1;
    }
    foreach my $gene (keys %keep_region) {
	push @{$hash{$gene}}, [$chrom,$start,$end];
    }
}

foreach $gname (keys %hash) {
    my $ct = 0;
    my @exons = @{$hash{$gname}};
    if ($strand{$gname} eq '-') {
	@exons = sort {$b->[2] <=> $a->[2]} @exons;
    }
    foreach $e (@exons) {
	$ct ++;
	print OUT join("\t",@{$e},join("|",$gname,$ct)),"\n";
    }
}
