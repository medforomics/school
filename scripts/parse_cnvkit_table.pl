#!/usr/bin/perl -w
#parse_cnvkit_table.pl

my @keep = ("AKT2","ATM","AURKA","BAP1","BCL2L1","BCL6","BIRC3","BRCA2","CCND1","CCNE1","CD79B","CDK8","CDKN2A","CDKN2B","CEBPA","EGFR","ERBB2","FGF10","FGF14","FGF19","FGF3","FGF4","FLT1","FLT3","FOXP1","GATA6","GNA13","GNAS","IKZF1","IL7R","IRS2","KDM4C","KLHL6","KMT2A","KRAS","LYN","MCL1","MITF","MYC","NOTCH2","PIK3CA","PRKAR1A","PRKDC","PTEN","RB1","RICTOR","RUNX1T1","SDHA","TFDP1","TP53","ZNF217");

my %keep = map {$_ => 1} @keep;


my @cvncalls = @ARGV;
foreach my $file (@ARGV) {
    open IN, "<$file" or die $!;
    my $out = $file;
    $out =~ s/call.cns/cnvcalls.txt/;
    open OUT, ">$out" or die $!;
    print OUT join("\t","CHROM","START","END","LEN","CN","DEPTH",
		   "Probes","Weight","GeneIDs"),"\n";
    my $header = <IN>;
    while (my $line = <IN>) {
	chomp($line);
	my ($chr,$start,$end,$geneids,$log2,$cn,$depth,
	    $probes,$weight) = split(/\t/,$line);
	my %genes;
	my @ids = split(/;|,/,$geneids);
	foreach my $gid (@ids) {
	    my ($key,$value) = split(/=/,$gid);
	    if ($key eq 'ensembl_gn' || $key eq 'identifier') {
		$genes{$value} = 1 if $keep{$value};
	    }
	}
	my $newgeneids = join(";", keys %genes);
	my $len = sprintf("%.1f",($end-$start)/1000000);
	next if ($cn == 2) || scalar(keys %genes) < 1;
	print OUT join("\t",$chr,$start,$end,$len,$cn,$depth,
		       $probes,$weight,$newgeneids),"\n";
    }
    close IN;
    close OUT;
}
