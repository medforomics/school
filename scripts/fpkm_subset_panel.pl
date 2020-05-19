#!/bin/bin/perl
use strict;
use warnings;
use diagnostics;
use Getopt::Long;


my ($genelist,$fpkm_file) = ("") x 2;

GetOptions( 'g|genes=s' => \$genelist,
            'f|fpkm=s' => \$fpkm_file);

my $fpkm_out = $fpkm_file;
$fpkm_out =~ s/\.fpkm.+txt/\.fpkm.capture.txt/;
open OUT, ">$fpkm_out" or die $!;

#gene_info.human.txt file is required for gene alias information
open GENEINFO, "</project/shared/bicf_workflow_ref/human/gene_info.human.txt" or die $!;
my %geneinfo;
while(my $gline = <GENEINFO>){
	chomp $gline;
	my ($taxid,$geneid,$symbol,$locus,$syn,@info) = split("\t",$gline);
	if($syn ne '-'){
		my @synonyms = split(/\|/,$syn);
		$geneinfo{$symbol} = $syn;
		foreach my $synonym(@synonyms){
			$geneinfo{$synonym} = $symbol;
		}
	}
}

#genelist is the list of genes in panel. Not specifying genelist will print out all lines(e.g. wholernaseq)
my %genes;
if($genelist ne ""){
	open GENES, "<$genelist" or die $!;
	while (my $gene = <GENES>){
		chomp $gene;
		$genes{$gene} = 1;
	}
	close GENES;
}
my $genelist_size = keys %genes;


#FPKM data
open FPKM, "<$fpkm_file" or die $!;
my %fpkmdata;
my %found;
my $header = <FPKM>;
chomp $header;
print OUT $header."\n";
while(my $line = <FPKM>){
	chomp $line;
	my ($gid,$gname,$ref,$strand,$start,$end,$cov,$fpkm,$tpm) = split("\t",$line);
	if($ref =~ /^chr/){ 
		if(!exists $fpkmdata{$gname}){$fpkmdata{$gname} = $line;}
		else{$fpkmdata{$gname} = $fpkmdata{$gname}."\n".$line;}	
	}
	if ($genelist_size == 0){ #wholernaseq
		print OUT $line."\n";
	}
	elsif($genes{$gname} and $ref =~ /^chr/){ #panel gene found in fpkm data
		print OUT $line."\n";
		$found{$gname} =1;
	}
}

foreach my $gene_key (keys %genes){
	if($found{$gene_key}){}
	else { #gene from gene list is not in fpkm data or has a different name
		if($geneinfo{$gene_key}){ #gene is in gene_info.human.txt
			if($fpkmdata{$geneinfo{$gene_key}}){ #gene from gene list is an alias in gene_info.human.txt
				print OUT $fpkmdata{$geneinfo{$gene_key}}."\n";
			}
			elsif($geneinfo{$gene_key} =~ /\|/){
				my @alias_array = split(/\|/, $geneinfo{$gene_key});
				foreach my $alias(@alias_array){
					if($fpkmdata{$alias}){print OUT $fpkmdata{$alias}."\n";}#gene from genelist is a Gene Name in gene_info.human.txt but fpkm data has an alias
				}
			}
		}
	}
}
close FPKM;
close OUT;
