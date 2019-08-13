#!/usr/bin/perl -w
#validation_json2txt.pl

use JSON;
my $jsonfile = shift @ARGV;
$jsontxt = `cat $jsonfile`;
$jsonref  = decode_json($jsontxt);

$prefix = (split(/\./,$jsonfile))[0];

my @vtypes = ('variants','cnvs','translocations');

foreach my $vtype (@vtypes) {
    my %csam;
    open OUT, ">$prefix\.$vtype\.txt" or die $!;
    foreach $tref (@{$jsonref->{$vtype}}) {
	%hash = %{$tref};
	next unless ($hash{selected});
	if ($vtype eq 'variants') {
	    print OUT join("\t",$prefix,$hash{'caseId'},$hash{chrom},$hash{pos},$hash{reference},$hash{alt},$hash{'tumorTotalDepth'},sprintf("%.2f",$hash{'tumorAltFrequency'}),
			   $hash{'normalTotalDepth'},sprintf("%.2f",$hash{'normalAltFrequency'}),$hash{'geneName'},$hash{'notation'}),"\n";
	}elsif ($vtype eq 'cnvs') {
	    my @genes = @{$hash{genes}};
	    foreach my $gid (@genes) {
		print OUT join("\t",$prefix,$hash{'caseId'},$gid,$hash{'copyNumber'},$hash{'aberrationType'},$hash{cytoband},$hash{score}),"\n";
	    }
	}else {
	    print OUT join("\t",$prefix,$hash{'caseId'},$hash{"fusionName"},$hash{"rnaReads"}),"\n";
	}
    }
    close OUT;
}

