#!/usr/bin/perl -w
#update_infofield_philips.pl

my $vcf = shift @ARGV;
open IN, "gunzip -c $vcf |" or die $!;
$outfile = $vcf;
$outfile =~ s/vcf.gz/test.vcf/;
open OUT, ">$outfile" or die $!;
while (my $line = <IN>) {
    chomp($line);
    if ($line =~ m/^#/) {
	print OUT $line,"\n";
    }else {
	my ($chrom, $pos,$id,$ref,$alt,$score,
	    $filter,$annot,$format,@gts) = split(/\t/, $line);
	my %hash = ();
	if ($pos eq 32337923) {
	    warn "Stop Here\n";
	}
	foreach $a (split(/;/,$annot)) {
	    my ($key,$val) = split(/=/,$a);
	    $hash{$key} = $val unless ($hash{$key});
	}
	my @deschead = split(/:/,$format);
	my $allele_info = shift @gts;
	@ainfo = split(/:/, $allele_info);
	foreach my $k (0..$#ainfo) {
	    $hash{$deschead[$k]} = $ainfo[$k];
	}
	my %fail;
	my $exacaf = '';
	if ($hash{AC_POPMAX} && $hash{AN_POPMAX}) {
	    @exacs = split(/,/,$hash{AC_POPMAX});
	    my $ac = 0;
	    foreach $val (@exacs) {
		$ac += $val if ($val =~ m/^\d+$/);
	    }
	    @exans = split(/,/,$hash{AN_POPMAX});
	    my $an = 0;
	    foreach $val (@exans) {
		$an += $val if ($val =~ m/^\d+$/);
	    }
	    $exacaf = sprintf("%.4f",$ac/$an) if ($ac > 0 && $an > 10);
	}
	unless ($exacaf eq '' || $exacaf <= 0.01) {
	    $fail{'COMMON'} = 1;
	}
	my $altct = 0;
	if ($hash{AO} =~ m/,/) { 
	    foreach (split(/,/, $hash{AO})) {
		$altct += $_;
	    }
	}else {
	    $altct = $hash{AO};
	}
	if ($hash{DP} =~ m/,/) {
	    @all = split(/,/, $hash{AD});
	    $hash{DP} = pop @all;
	}
	$hash{AF} = sprintf("%.3f",$altct/$hash{DP});
	my $cosmicsubj = 0;
	@cosmicct = split(/,/,$hash{CNT}) if $hash{CNT};
	foreach $val (@cosmicct) {
	    $cosmicsubj += $val if ($val =~ m/^\d+$/);
	}
	if ($hash{DP} < 20) {
	    $fail{'LowDepth'} = 1;
	}
	if (($hash{CallSet} =~ m/hotspot/ || $id =~ m/COS/) && $cosmicsubj >= 5) {
	    $fail{'LowAltCt'} = 1 if ($altct < 3);
	    $fail{'LowMAF'} = 1 if ($hash{AF} < 0.01);
	}else {
	    $fail{'OneCaller'} = 1 if ($hash{CallSet} =~ m/,/);
	    $fail{'LowAltCt'} = 1 if ($altct < 8);
	    $fail{'LowMAF'} = 1 if ($hash{AF} < 0.05);
	}
	
	my @nannot;
	foreach $info (sort {$a cmp $b} keys %hash) {
	    if ($hash{$info}) {
		push @nannot, $info."=".$hash{$info} 
	    }else {
		push @nannot, $info;
	    }
	}
	my @fail = keys %fail;
	if (scalar(@fail) == 0) {
	    $filter = 'PASS';
	}else {
	    $filter = join(";", 'FailedQC',@fail);
	}
	$newannot = join(";",@nannot);
	print OUT join("\t",$chrom, $pos,$id,$ref,$alt,$score,
		       $filter,$newannot,$format,$allele_info),"\n";
    }
}
close OUT;
close IN;
