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
	my ($ensid,$gene,$chr,$strand,$start,$end,$cov,$fpkm,$tmp) = split(/\t/,$line);
	$cts{$gene}{$sample} = $fpkm;
    }
    close IN;
}

my %exc = (HBB=>1,HBD=>1,HBBP1=>1,
	   HBG1=>1,HBG2=>1,HBE1=>1,
	   HBZ=>1,HBZP1=>1,HBA2=>1,HBA1=>1);

open OUT, ">$opt{outdir}\/countTable.fpkm.txt" or die $!;
my @samples = sort {$a cmp $b} keys %samps;
print OUT join("\t",'ENSEMBL','SYMBOL','TYPE',@samples),"\n";

foreach $gene (keys %cts) {
    my @line;
    foreach $sample (@samples) {
	unless ($cts{$gene}{$sample}) {
	    $cts{$gene}{$sample} = 0;
	}
	push @line, $cts{$gene}{$sample};
    }
    my @sortct = sort {$b cmp $a} @line;
    next unless ($sortct[0] >= 1);
    next if ($names{$gene}{type} eq 'rRNA');
    next if ($exc{$gene});
    print OUT join("\t",$names{$gene}{ensembl},$gene,
		   $names{$gene}{type},@line),"\n";
}
close OUT;
