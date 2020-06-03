#!/usr/bin/perl 
#migrate_db.pl

my ($caseid,$tumorid,$dnarunid) = @ARGV;

my $prefix = "$caseid\/$caseid";
my $passvcf = "$caseid\/$caseid\.pass.vcf.gz";

my %info;
my %ref;
open PASS, "gunzip -c $passvcf |" or die $!;
while (my $line = <PASS>) {
    chomp($line);
    next if ($line =~ m/^#/);
    my ($chrom,$pos,$id,$ref,$alt,$score,$filter,$annot,@info) = split(/\t/,$line);
    my %hash = ();
    foreach $a (split(/;/,$annot)) {
	my ($key,$val) = split(/=/,$a);
	$val =~ s/,/\|/g if ($val);
	$hash{$key} = $val unless ($hash{$key});
    }
    my $key = join(":",$chrom,$pos,$ref,$alt);
    my $indellen = abs(length($ref) - length($alt));
    $ref{$chrom}{$pos}{$ref}{$alt} = $key;
    $info{$key} = {tumoraf=>$hash{AF},tumordp=>$hash{DP},
		   tumorao=>sprintf("%.0f",$hash{DP}*$hash{AF}),
		   normalaf=>$hash{NormalAF},normaldp=>$hash{NormalDP},
		   normalao=>sprintf("%.0f",$hash{NormalDP}*$hash{NormalAF}),
		   indellen=>$indellen};
    $info{$key}{FS} = $hash{FS} if ($hash{FS});
}
my %winbed;
open WINBED, "<$prefix\.50window.bed" or die $!;
while (my $line = <WINBED>) {
    chomp($line);
    my ($chrom,$start,$end,$variant) = split(/\t/,$line);
    $gckey = $chrom.":".$start."-".$end;
    $win50{$gckey} = $variant
}
open GC, "<$prefix\.50window.gc.txt" or die $!;
while (my $line = <GC>) {
    chomp($line);
    my ($key,$seq,$a,$gc) = split(/\t/,$line);
    my $variant = $win50{$key};
    $info{$variant}{gcontent} = $gc;
}
close GC;
my %hp;
open HP, "<$prefix\.50window.homopolymer.txt" or die $!;
while (my $line = <HP>) {
    chomp($line);
    my ($key,$pname,$pattern,$strand,$start,$end,$match) = split(/\t/,$line);
    foreach ($start..$end) {
	$hp{$key}{$_} = 1;
    }
}
close HP;
foreach my $reg (keys %hp) {
    my $variant = $win50{$reg};
    my ($chr,$start,$end) = split(/:|-/,$reg);
    my $hpnum = scalar(keys %{$hp{$reg}});
    $info{$variant}{homopolymer} = $hpnum;
}
open RPEAT, "<$prefix\.repeat.vcf.txt" or die $!;
while (my $line = <RPEAT>) {
    chomp($line);
    my ($chr,$start,$end,$key,$chr2,$s1,$e2,$rtype,$len) = split(/\t/,$line);
    my $varlen = $end-$start;
    $percrep = sprintf("%.2f",$len/$varlen);
    $info{$key}{repeatperc} = $percrep;
}
close RPEAT;
open BED, "<$prefix\.bed" or die $!;
while (my $line = <BED>) {
    chomp($line);
    my ($chrom,$start,$end,$variant) = split(/\t/,$line);
    $statkey = $chrom.":".$start."-".$end;
    $winbed{$statkey} = $variant
}

open BSTAT, "<$prefix\.bamstat" or die $!;
my $header = <BSTAT>;
chomp($header);
my @colnames = split(/\t/,$header);
while (my $line = <BSTAT>) {
    chomp($line);
    my @row = split(/\t/,$line);
    my %hash;
    foreach my $i (0..$#colnames) {
 	$hash{$colnames[$i]} = $row[$i];
    }
    my $statkey = $hash{chrom}.":".$hash{chromStart}."-".$hash{chromEnd};
    my $variant = $winbed{$statkey};
    foreach my $metric (keys %hash) {
	next if ($metric =~ m/chrom/);
	$info{$variant}{$metric} = $hash{$metric};
    }
}
close BSTAT;
open BAMCT, "<$prefix\.bamreadct.txt" or die $!;
my %bq;
while(my $line=<BAMCT>){
    chomp($line);
    my ($chrom,$pos,$ref,$depth,@counts) = split("\t",$line);
    foreach my $ct(@counts){
	my ($alt,$depth,$mapqual,$basequal,@ostat) = split(/:/,$ct);
	if ($ref{$chrom}{$pos}{$ref}{$alt}) {
	    $key = join(":",$chrom,$pos,$ref,$alt);
	    $info{$key}{basequal}=$basequal;
	}
    }
}
my @statfile = `ls $caseid\/dna*\/*.stat.txt $caseid\/somatic*\/*.stat.txt`;
foreach my $txt (@statfile) {
    open SFILE, "<$txt" or die $!;
    my $header = <SFILE>;
    chomp($header);
    my @colnames = split(/\t/,$header);
    while (my $line = <SFILE>) {
	chomp($line);
	my @row = split(/\t/,$line);
	my %hash;
	foreach my $i (0..$#colnames) {
	    $hash{$colnames[$i]} = $row[$i];
	}
	my $key = join(":",$hash{CHROM},$hash{POS},$hash{REF},$hash{ALT});
	next unless $info{$key};
	$info{$key}{MapQual} = $hash{MQ};
	$info{$key}{MapPOS} = (split(/,/,$hash{MPOS}))[1];
	$info{$key}{MapQualAlt} = (split(/,/,$hash{MMQ}))[1];
	$info{$key}{MapQualRef} = (split(/,/,$hash{MMQ}))[0];
	$info{$key}{FS} = $hash{FS};
    }
} 
open OUT, ">$prefix\.vcfstats.txt" or die $!;
F1:foreach $key (sort {$a cmp $b} keys %info) {
    my @nannot;
    foreach $statid (sort {$a cmp $b} keys %{$info{$key}}) {
	if (exists $info{$key}{$statid}) {
	    push @nannot, $statid."=".$info{$key}{$statid};
	}
    }
    push @nannot, "CaseID=".$caseid;
    my $newannot = join(";",@nannot);
    next unless $key;
    print OUT join("\t",split(/:/,$key),$newannot),"\n";
}
close OUT;
