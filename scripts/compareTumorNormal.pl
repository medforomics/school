#!/usr/bin/perl -w
#integrate_datasets.pl

my @files = @ARGV;

foreach my $file (@files) {
    open IN, "gunzip -c $file |" or die $!;
    my $total = 0;
    my $rnatotal = 0;
    my %concordTumor;
    my $concordRNASeq = 0;
    
    while (my $line = <IN>) {
	chomp($line);
	next if ($line =~ m/#/);
	my ($chr, $pos,$id,$ref,$alt,$score,
	    $filter,$annot,$format,@gts) = split(/\t/, $line);
	next unless ($filter eq 'PASS' || $filter eq 'FailedQC;COMMON');
	my %hash = ();
	foreach $a (split(/;/,$annot)) {
	    my ($key,$val) = split(/=/,$a);
	    $hash{$key} = $val unless ($hash{$key});
	}
	my @deschead = split(/:/,$format);
	my $allele_info = shift @gts;
	@ainfo = split(/:/, $allele_info);
	my %gtinfo = ();
	my @mutallfreq = ();
	foreach my $k (0..$#ainfo) {
	    $gtinfo{$deschead[$k]} = $ainfo[$k];
	}
	next unless ($hash{SS} == 1);
	$total ++;
	my $normalAF = (split(/,/,$hash{NormalAF}))[0];
	my $tumorAF = (split(/,/,$hash{AF}))[0];
	my @diffs = (0.2,0.3);
	my @tumornums;
	foreach $diffthresh (@diffs) {
	    $sampdiff = abs($normalAF-$tumorAF);
	    if ($tumorAF  && $sampdiff <= $diffthresh) { 
		$concordTumor{$diffthresh} ++;
	    }
	}
	if ($hash{RnaSeqDP} && $hash{RnaSeqDP} > 20) {
	    $rnatotal ++;
	    $hash{RnaSeqAF} = 0 unless ($hash{RnaSeqAF});
	    $rnaseqAF = (split(/,/,$hash{RnaSeqAF}))[0];
	    if ($rnaseqAF > 0.05) { 
		$concordRNASeq ++;
	    }else {
		#warn "discordant\n";
	    }
	}
    }
    my $perc_concordRNASeq = 0;
    $perc_concordRNASeq = sprintf("%.3f",$concordRNASeq/$rnatotal) if ($rnatotal > 1);
    print join("\t",$file,$total,sprintf("%.3f",$concordTumor{0.3}/$total),
	       $rnatotal,$perc_concordRNASeq),"\n";
}
