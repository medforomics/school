use strict;
use warnings;
use List::Util qw(min max);

my ($inputFile, $geneInfoFile, $refversion, @databaseList) = @ARGV;
my @databaseOptionList = map {"--$_"} @databaseList;
my %symbolHash = ();
my %synonymSymbolHash = ();
{
	open(my $reader, $geneInfoFile);
	my @columnList = ();
	while(my $line = <$reader>) {
		chomp($line);
		if($line =~ s/^#//) {
			@columnList = split(/\t/, $line);
		} else {
			my %tokenHash = ();
			@tokenHash{@columnList} = split(/\t/, $line);
			$symbolHash{$tokenHash{'Symbol'}} = 1;
			if($tokenHash{'Synonyms'} ne '-') {
				$synonymSymbolHash{$_} = $tokenHash{'Symbol'} foreach(split(/\|/, $tokenHash{'Synonyms'}));
			}
		}
	}
	close($reader);
}
my %geneRegionListHash = ();
my %geneRegionListListHash = ();
{
	chomp(my @configList = `transvar current --refversion $refversion`);
	my %configHash = map {$_->[0] => $_->[1]} map {[split(/: /, $_, 2)]} @configList;
	open(my $reader, "cat @configHash{@databaseList} |");
	while(my $line = <$reader>) {
		chomp($line);
		my @tokenList = split(/\t/, $line);
		push(@{$geneRegionListHash{$tokenList[0]}}, "$tokenList[6]:$tokenList[4]-$tokenList[5]");
		if($tokenList[10] =~ /^\[\((.+)\)\]$/) {
			my $regions = $1;
			my @regionList = map {/^([0-9]+), ([0-9]+)$/ ? "$tokenList[6]:$1-$2" : ()} split(/\), \(/, $regions);
			@regionList = reverse @regionList if($tokenList[7] eq '-');
			push(@{$geneRegionListListHash{$tokenList[0]}}, \@regionList);
		}
	}
	close($reader);
}
open(my $reader, $inputFile);
while(my $line = <$reader>) {
	chomp($line);
	s/^ *//, s/ *$// foreach(my ($sample, $gene, $mutation) = split(/\t/, $line));
	if(defined($symbolHash{$gene})) {
	} elsif(defined($gene = $synonymSymbolHash{$gene})) {
	} else {
		print STDERR join("\t", 'unknown_gene', $line), "\n";
		next;
	}
	if($mutation =~ /^[A-Z][0-9]+[A-Z]$/) { # missense
		printBedLine($line, getRegionList(transvar('p', $gene, $mutation)));
	} elsif($mutation =~ /^[A-Z][0-9]+$/) { # aminoacid
		printBedLine($line, getRegionList(transvar('p', $gene, $mutation)));
	} elsif($mutation =~ /^[A-Z][0-9]+[*]$/) { # nonsense
		printBedLine($line, getRegionList(transvar('p', $gene, $mutation)));
	} elsif($mutation =~ /^[*][0-9]+[A-Z]$/) { # readthrough
		printBedLine($line, getRegionList(transvar('p', $gene, $mutation)));
	} elsif($mutation =~ /^[A-Z][0-9]+_[A-Z][0-9]+/) { # deletion or insertion
		printBedLine($line, getRegionList(transvar('p', $gene, $mutation)));
	} elsif($mutation =~ /^[A-Z][0-9]+del/) { # aminoacid deletion
		printBedLine($line, getRegionList(transvar('p', $gene, $mutation)));
	} elsif($mutation =~ /^[A-Z][0-9]+fs/) { # frameshift
		printBedLine($line, getRegionList(transvar('p', $gene, $mutation)));
	} elsif($mutation =~ /^[0-9]+[+-][0-9]+/) { # splicing
		printBedLine($line, getRegionList(transvar('c', $gene, $mutation)));
	} elsif($mutation eq 'amplification' || $mutation eq 'rearrangement') { # amplification or rearrangement
		printBedLine($line, defined($_ = $geneRegionListHash{$gene}) ? @$_ : ());
	} elsif($mutation =~ /^loss exons ([0-9]+)-([0-9]+)$/) { # loss exons
		my @regionList = ();
		foreach(map {[$_->[$1 - 1], $_->[$2 - 1]]} @{$geneRegionListListHash{$gene}}) {
			my ($region1, $region2) = @$_;
			if(defined($region1) && defined($region2)) {
				my ($chromosome1, $start1, $end1) = ($region1 =~ /^(.+):([0-9]+)-([0-9]+)$/) ? ($1, $2, $3) : ('', '', '');
				my ($chromosome2, $start2, $end2) = ($region2 =~ /^(.+):([0-9]+)-([0-9]+)$/) ? ($1, $2, $3) : ('', '', '');
				if($chromosome1 eq $chromosome2) {
					my $region = sprintf('%s:%d-%d', $chromosome1, min($start1, $start2), max($end1, $end2));
					push(@regionList, $region);
				}
			}
		}
		printBedLine($line, @regionList);
	} else {
		print STDERR join("\t", 'unknown_mutation', $line), "\n";
	}
}
close($reader);

sub printBedLine {
	my ($line, @regionList) = @_;
	if(@regionList) {
		open(my $writer, "| sort --field-separator='\t' -k1,1 -k2,2n -k3,3n | uniq");
		print $writer join("\t", @$_, $line), "\n" foreach(map {/^(.+):([0-9]+)-([0-9]+)$/ ? [$1, $2 - 1, $3] : ()} @regionList);
		close($writer);
	} else {
		print STDERR join("\t", 'no_region', $line), "\n";
	}
}

sub getRegionList {
	my (@tokenHashList) = @_;
	my @regionList = ();
	foreach(@tokenHashList) {
		my %tokenHash = %$_;
		if(defined($_ = $tokenHash{'coordinates(gDNA/cDNA/protein)'}) && (my @coordinateList = split(/\//, $_))) {
			if($coordinateList[0] =~ /^(.+):g.([0-9]+)_([0-9]+)/) {
				push(@regionList, "$1:$2-$3");
			} elsif($coordinateList[0] =~ /^(.+):g.([0-9]+)/) {
				push(@regionList, "$1:$2-$2");
			}
		}
	}
	return @regionList;
}

sub transvar {
	my ($type, $gene, $mutation) = @_;
	my @tokenHashList = ();
	open(my $reader, "transvar ${type}anno --refversion $refversion -i '$gene:$type.$mutation' @databaseOptionList 2> /dev/null |");
	my @columnList = ();
	{
		chomp(my $line = <$reader>);
		@columnList = split(/\t/, $line);
	}
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		@tokenHash{@columnList} = split(/\t/, $line);
		push(@tokenHashList, \%tokenHash);
	}
	close($reader);
	return @tokenHashList;
}
