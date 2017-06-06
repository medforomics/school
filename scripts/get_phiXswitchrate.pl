#!/usr/bin/perl -w
#get_phiXswitchrate.pl

open CTS, "<phiXcounts.txt" or die $!;
while (my $line = <CTS>) {
    chomp($line);
    my ($sampleName,$phiXreads) = split(/\t/,$line);
    $SampleMergeName = $sampleName;
    my @samplename = split(/_/,$sampleName);
    if ($samplename[-1] =~ m/^[A|B|C|D]$/ || $samplename[-1] =~ m/^[A|T|C|G]+$/) {
	pop @samplename;
	$SampleMergeName = join("_",@samplename);
    }
    $ct{$SampleMergeName} += $phiXreads;
}
my %info;
open SNAME, "<sampnamebyflowcell.txt" or die $!;
while (my $line = <SNAME>) {
    chomp($line);
    my ($flowcell, $samplename,$demultiplexname,$lane,$clusters) = split(/\t/,$line);
    next if ($samplename =~ m/exome/i || $flowcell =~ m/M70356|160930/);
    unless ($ct{$samplename}) {
	print $line,"\n";
    }else {
	$info{$flowcell}{total} += $clusters;
	$info{$flowcell}{phiX} += $ct{$samplename};
    }
}

open SRUN, "<seqrun.txt" or die $!;
open OUT, ">seqrun.info.txt" or die $!;
while (my $line = <SRUN>) {
    chomp($line);
    my @ary = split(/\t/,$line);
    if ($info{$ary[0]}) {
	print OUT join("\t",$line,sprintf("%.1e",$info{$ary[0]}{phiX}/$info{$ary[0]}{total})),"\n";
    }else {
	#print $line,"\n";
    }
}
