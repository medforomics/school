#!/usr/bin/perl -w
#archer2answer.pl

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my %opt = ();
my $results = GetOptions (\%opt,'fusion|f=s','prefix|p=s','help|h');
open OAS, ">$opt{prefix}\.translocations.answer.txt" or die $!;
print OAS join("\t","FusionName","LeftGene","LefttBreakpoint","LeftGeneExons","LeftStrand",
	       "RightGene","RightBreakpoint","RightGeneExons","RightStrand",
	       "RNAReads","DNAReads","FusionType","Annot",'Filter','ChrType','ChrDistance','PercentSupportingReads'),"\n";

open FUSION, "<$opt{fusion}" or die $!;
my $header = <FUSION>;
chomp($header);
$header =~ s/\"//g;
my @hline = split(/,/,$header);
while (my $line = <FUSION>) {
    chomp($line);
    my @row = split(/,/,$line);
    my %hash;
    foreach my $i (0..$#row) {
	$hash{$hline[$i]} = $row[$i];
    }
    my @annot2;
    my @filter;
    if ($hash{'KNOWN QUIVER FUSION'} ne 'FALSE') {
	push @annot2, 'QUIVER';
    }if ($hash{'LOW CONFIDENCE FUSION'} ne 'FALSE') {
	push @filter, 'LowConfidenceFusion';
    }if ($hash{'KNOWN ENSEMBLE PARALOG'} ne 'FALSE') {
	push @filter, 'KnownEnsembleParalog';
    }if ($hash{'TRANSCRIPT READTHROUGH'} ne 'FALSE') {
	push @filter, 'ReadThrough';
    }if ($hash{'MISPRIMING EVENT'} ne 'FALSE') {
	push @filter, 'MisprimingEvent';
    }if ($hash{'TOTAL NUMBER OF SUPPORTIVE READS'} < 50) {
	push @filter, 'LowReadCt';
    }if ($hash{'NUMBER OF SUPPORTIVE READS WITH UNIQUE START SITES'} < 2) {
	push @filter, 'LowUniqueStartReadCt';
    }if ($hash{"PERCENT OF READS AT BREAKPOINT SUPPORTING FUSION"} < 10) {
	push @filter, 'LowPercentReads';
    }
    my ($left_chr,$left_pos,$left_strand) = split(/:/,$hash{"GENE 1 BREAKPOINT"});
    $hash{LeftBreakpoint} = join(":",$left_chr,$left_pos);
    $hash{LeftStrand} = $left_strand;
    my ($leftgeneids,$lexon,$lloci) = split(/\|/,$hash{"5PPRIME EXON INVOLVED IN FUSION"});
    $lexon =~ m/exon:(\d+)/;
    $leftexon = $1; 
    my ($leftgene,$leftgenestrand,$leftgene_acc) = split(/\(|\)/,$leftgeneids);
    $hash{LeftGene} = $leftgene;
    my ($right_chr,$right_pos,$right_strand) = split(/:/,$hash{"GENE 1 BREAKPOINT"});
    $hash{RightBreakpoint} = join(":",$right_chr,$right_pos);
    $hash{RightStrand} = $right_strand;
    my ($rightgeneids,$rexon,$rloci) = split(/\|/,$hash{"5PPRIME EXON INVOLVED IN FUSION"});
    $rexon =~ m/exon:(\d+)/;
    $rightexon = $1; 
    my ($rightgene,$rightgenestrand,$rightgene_acc) = split(/\(|\)/,$rightgeneids);
    $hash{RightGene} = $rightgene;

    $fusion_annot = '';
    if (scalar(@annot2) > 0) {
	$fusion_annot = join(",",@annot2);
    }
    my $qc ='PASS';
    if (scalar(@filter) > 0) {
	$qc = join(";","FailedQC",@filter);
    }
    if ($left_chr eq $right_chr) {
	$chrtype = 'INTRACHROMOSOMAL';
	$diff = sprintf("%.2f",abs($right_pos-$left_pos)/1000000);
	$chrdist = join(":",$left_chr,$diff."Mb");
    }else {
	$chrtype = 'INTERCHROMOSOMAL';
	$chrdist = join("--",$left_chr,$right_chr);
    }

    print OAS join("\t",$hash{FUSION},$hash{LeftGene},$hash{LeftBreakpoint},$leftexon,$hash{LeftStrand},
		   $hash{RightGene},$hash{RightBreakpoint},$rightexon,$hash{RightStrand},
		   $hash{"TOTAL NUMBER OF SUPPORTIVE READS"},$hash{"NUMBER OF SUPPORTIVE READS WITH UNIQUE START SITES"},'',$fusion_annot,$qc,$chrtype,$chrdist,$hash{"PERCENT OF READS AT BREAKPOINT SUPPORTING FUSION"}),"\n";
}
