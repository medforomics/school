#!/usr/bin/perl -w
#parse_viral_idxstats.pl

my %genoname;
open GNAME, "<viral_genome_names.txt" or die $!;
while (my $line = <GNAME>) {
    chomp($line);
    my ($short_name,$long_name) = split(/\t/,$line);
    my @longname = split(/_/,$long_name);
    if ($long_name =~ m/^NC_/) {
	$genoname{join("_",@longname[0..1])} = $short_name;
    }else {
	$genoname{$longname[0]} = $short_name;
    }
}

open OUT, ">viral_results.txt" or die $!;

my @idxstats = @ARGV;
foreach $file (@idxstats) {
    open IN, "<$file" or die $!;
    my ($case,$sample,$filename) = split(/\//,$file);
    while (my $line = <IN>) {
	chomp($line);
	my ($geno,$length,$mapped,$unmapped) = split(/\t/,$line);
	$geno = (split(/\./,$geno))[0];
	print OUT join("\t",$case,$sample,$genoname{$geno},$geno,$mapped),"\n" if ($mapped > 10);
    }
}
