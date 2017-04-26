#!/usr/bin/perl

use strict;

if (scalar @ARGV != 1) {
	print "vcf2bed.pl <vcf_file>\n";
	exit 1;
}

my $vcfFile = $ARGV[0];

open VCF, $vcfFile
	or die "Cannot open $vcfFile";

my $ct = 0;
while (<VCF>) {
    next if (index($_, '#') == 0);
    my @fields = split /\t/, $_;
    $ct ++;
    my $chrom = $fields[0];
    my $pos = $fields[1];
    my $ref = $fields[3];
    my $annot = $fields[7];
    my $alt = $fields[4];
    if ($fields[2] eq 'N') {
	$fields[2] = 'NB'.sprintf("%06s",$ct);
    }
    my %hash = ();
    next if ($chrom =~ m/_/); 
    foreach $a (split(/;/,$annot)) {
	my ($key,$val) = split(/=/,$a);
	$hash{$key} = $val unless ($hash{$key});
    }
    $hash{'END'} = $pos+1 unless $hash{'END'};
    next if ($hash{CHR2} && $hash{CHR2} =~ m/_/); 
    if ($alt =~ m/chr(\w+):(\d+)/i) {
	if ($1 eq $chrom) {
	    print join("\t",$chrom,$pos,$2,$fields[2]),"\n";
	}elsif ($fields[2] =~ m/_\d+/) {
	    print join("\t",$chrom,$pos,$hash{END},$fields[2]),"\n";
	}else {
	    print join("\t",$chrom,$pos,$hash{END},$fields[2]."_1"),"\n";
	    print join("\t",'chr'.$1,$2,$2+1,$fields[2]."_2"),"\n";
	}
    }else {
	if ($hash{CHR2} && $hash{CHR2} eq $chrom) {
	    print join("\t",$chrom,$pos,$hash{END},$fields[2]),"\n";
	}elsif ($hash{CHR2} && $hash{CHR2} ne $chrom) {
	    print join("\t",$chrom,$pos,$pos+1,$fields[2]."_1"),"\n";
	    print join("\t",$hash{CHR2},$hash{END},$hash{END}+1,$fields[2]."_2"),"\n";
	}unless ($hash{CHR2}) {
	    print join("\t",$chrom,$pos,$hash{END},$fields[2]),"\n";
	}
    }
}

close VCF;
