#!/usr/bin/perl -w
#concat_edgeR.pl

#symbol	chrom	start	end	ensembl	type	logFC	logCPM	PValue	monocytes	neutrophils	rawP	fdr	bonf

my @files = @ARGV;
open OUT, ">edgeR.results.txt" or die $!;
print OUT join("\t",'G1','G2','ENSEMBL','SYMBOL','TYPE','logFC','logCPM',
	       'pval','fdr','G1.Mean','G2.Mean'),"\n";
foreach $f (@files) {
    open IN, "<$f" or die $!;
    my ($fname,$ext) = split(/\.edgeR\./,$f);
    $fname =~ s/edgeR\///;
    my ($g1,$g2) = split(/_/,$fname);
    my $head = <IN>;
    chomp($head);
    my @colnames = split(/\t/,$head);
    while (my $line = <IN>) {
	chomp($line);
	my @row = split(/\t/,$line);
	my %hash;
	foreach my $i (0..$#row) {
	    $hash{$colnames[$i]} = $row[$i];
	}
	print OUT join("\t",$g1,$g2,$hash{ensembl},$hash{symbol},$hash{type},
		      sprintf("%.2f",$hash{logFC}),
		       sprintf("%.2f",$hash{logCPM}),
		       sprintf("%.1e",$hash{rawP}),sprintf("%.1e",$hash{fdr}),
		       sprintf("%.0f",$hash{$g1}),sprintf("%.0f",$hash{$g2})),"\n";
    }
}
