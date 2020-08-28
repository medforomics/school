#!/usr/bin/perl -w
#parse_viral_idxstats.pl

my %genoname;
open GNAME, "</project/shared/bicf_workflow_ref/human_virus_genome/clinlab_idt_genomes/viral_genome_names.txt" or die $!;
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
print OUT join("\t","SampleID","VirusName","VirusAcc","VirusDescription","ViralReadCt"),"\n";

my @idxstats = @ARGV;
foreach $file (@idxstats) {
    open IN, "<$file" or die $!;
    my ($dir,$filename) = split(/\//,$file);
    my ($sample,@dir) = split(/\./,$filename);
    my $header = <IN>;
    chomp($header);
    my @colnames = split(/\t/,$header);
    while (my $line = <IN>) {
	chomp($line);
	my @row = split(/\t/,$line);
	my %hash;
	foreach my $i (0..$#row) {
	    $hash{$colnames[$i]} = $row[$i];
	}
	$hash{VirusAcc} = (split(/\./,$hash{VirusAcc}))[0];
	print OUT join("\t",$sample,$genoname{$hash{VirusAcc}},$hash{VirusAcc},$hash{VirusName},$hash{Mapped}),"\n" if ($hash{Mapped} > 5);
    }
}
