#!/usr/bin/perl -w
#temp.pl

my $ss = shift @ARGV;
open SS, "<$ss" or die $!;
my $header = <SS>;
chomp($header);
my @colnames = split(/\t/,$header);

while (my $line = <SS>) {
    chomp($line);
    my @row = split(/\t/,$line);
    my %hash;
    foreach $i (0..$#row) {
		$hash{$colnames[$i]} = $row[$i];
    }
    next unless ($hash{Class} =~ m/Tumor|Normal/i);
    
    $samps{$hash{SubjectID}}{lc($hash{Class})} = $hash{SampleMergeName};
}

open TONLY, ">design_tumor_only.txt" or die $!;
open TNPAIR, ">design_tumor_normal.txt" or die $!;
print TNPAIR join("\t",'TumorID','NormalID','TumorBAM','NormalBAM',
		  'TumorOntargetBAM','NormalOntargetBAM'),"\n";
print TONLY join("\t",'SampleID','BAM','OntargetBAM'),"\n";
foreach my $subjid (keys %samps) {
    my @ctypes = keys %{$samps{$subjid}};
    if (scalar(@ctypes) > 1) {
	print TNPAIR join("\t",$samps{$subjid}{tumor},$samps{$subjid}{normal},
			  $samps{$subjid}{tumor}.".bam",$samps{$subjid}{normal}.".bam",
			  $samps{$subjid}{tumor}.".ontarget.bam",
			  $samps{$subjid}{normal}.".ontarget.bam"),"\n";
    } else {
	print TONLY join("\t",$samps{$subjid}{$ctypes[0]},$samps{$subjid}{$ctypes[0]}.".bam",
			 $samps{$subjid}{$ctypes[0]}.".ontarget.bam"),"\n";
    }
}
close TNPAIR;
close TONLY;
