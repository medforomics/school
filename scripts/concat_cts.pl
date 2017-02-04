#!/usr/bin/perl -w
#merge_tables.pl

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use File::Basename;

my %opt = ();
my $results = GetOptions (\%opt,'genenames|g=s','outdir|o=s','help|h');

#$gnames = $opt{genenames};

open SYM, "<genenames.txt" or die $!;
my %symbs = ();
while (my $line = <SYM>) {
    chomp($line);
    my ($chrom,$start,$end,$ensembl,$symbol,$type) = split(/\t/,$line);
    $ensembl = (split(/\./,$ensembl))[0];
    $names{$symbol} = {ensembl=>$ensembl,type=>$type};
}

my @files = @ARGV;
foreach $file (@files) {
    chomp($file);
    open IN, "<$file" or die $!;
    $fname = basename($file);
    my $sample = (split(/\./,$fname))[0];
    $samps{$sample} = 1;
    my $command = <IN>;
    my $head = <IN>;
    chomp($head);
    my @colnames = split(/\t/,$head);
    while (my $line = <IN>) {
	chomp($line);
	my @row = split(/\t/,$line);
	my $gene = $row[0];
	my $ct = $row[-1];
	my $type = $names{$gene}{'type'};
	$type = $gene unless ($names{$gene}{'type'});
	$types{$type} = 1;
	$total{$sample} += $ct;
	$readcts{$sample}{$type} += $ct;
	next if($gene =~ m/^__/);
	$cts{$gene}{$sample} = $ct;
    }
    close IN;
}

my @gtypes = sort {$a cmp $b} keys %types;
open OUT, ">$opt{outdir}\/countTable.stats.txt" or die $!;
print OUT join("\t","Sample","Type","ReadPerc","ReadCt","TotalReads"),"\n";

foreach $sample (sort {$a cmp $b} keys %readcts) {
    my @line;
    my $total = 0;
    foreach $type (@gtypes) {
	$ct = 0;
	$ct = sprintf("%.2e",$readcts{$sample}{$type}/$total{$sample}) if ($readcts{$sample}{$type});
	print OUT join("\t",$sample,$type,$ct,$readcts{$sample}{$type},$total{$sample}),"\n";
    }
}

my %exc = (HBB=>1,HBD=>1,HBBP1=>1,
	   HBG1=>1,HBG2=>1,HBE1=>1,
	   HBZ=>1,HBZP1=>1,HBA2=>1,HBA1=>1);

open OUT, ">$opt{outdir}\/countTable.txt" or die $!;
open CPM, ">$opt{outdir}\/countTable.logCPM.txt" or die $!;
my @samples = sort {$a cmp $b} keys %samps;
print OUT join("\t",'ENSEMBL','SYMBOL','TYPE',@samples),"\n";
print CPM join("\t",'ENSEMBL','SYMBOL','TYPE',@samples),"\n";

foreach $gene (keys %cts) {
    my @line;
    my @cpm;
    foreach $sample (@samples) {
	unless ($cts{$gene}{$sample}) {
	    $cts{$gene}{$sample} = 0;
	}
	$cpm = ($cts{$gene}{$sample}/$total{$sample})*1e6;
	push @cpm, sprintf("%.2f",log2($cpm));
	push @line, $cts{$gene}{$sample};
    }
    my @sortct = sort {$b cmp $a} @cpm;
    next unless ($sortct[0] >= 1);
    next if ($names{$gene}{type} eq 'rRNA');
    next if ($exc{$gene});
    print OUT join("\t",$names{$gene}{ensembl},$gene,
		   $names{$gene}{type},@line),"\n";
    print CPM join("\t",$names{$gene}{ensembl},$gene,
		   $names{$gene}{type},@cpm),"\n";
}
close OUT;
close CPM;
sub log2 {
    $n = shift @_;
    if ($n < 1) {
	return 0;
    }else {
	return(log($n)/log(2));
    }
}
