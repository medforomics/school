#!/usr/bin/perl -w
#parse_getsampleJson.pl

use JSON;
my $jsonfile = shift @ARGV;
$jsontxt = `cat $jsonfile`;
$jsonref  = decode_json($jsontxt);

die unless $jsonref->{'limsId'};
$limsid = $jsonref->{'limsId'};
$projectid = $jsonref->{'projectId'};

my @normals;
my @tumors;
my @rnaseq;

my @stypes = ('nDnaSamples','tDnaSamples','tRnaSamples');
my %caseSamples;

foreach my $stype (@stypes) {
    my %csam;
    foreach $tref (@{$jsonref->{$stype}}) {
	%hash = %{$tref};
	next unless ($hash{runId});
	push @{$csam{$hash{analysisSample}}}, [$hash{runId},$hash{sampleId}];
    }
    if ($csam{1}) {
	$csam{'true'} = $csam{1};
    }if ($csam{0}) {
	$csam{'false'} = $csam{0};
    }
    if ($csam{'true'}) {
	my @samples = @{$csam{'true'}};
	die "too many ".$stype if (scalar(@samples) > 1);
	$caseSamples{$stype} = $samples[0];
    }elsif($csam{'false'}) {
	my @samples = @{$csam{'false'}};
	die "too many ".$stype if (scalar(@samples) > 1);
	$caseSamples{$stype} = $samples[0];
    }else {
	$caseSamples{$stype} = ['NA','NA'];
    }
}

print join(" ",$limsid,$projectid,@{$caseSamples{'tDnaSamples'}},@{$caseSamples{'nDnaSamples'}},@{$caseSamples{'tRnaSamples'}}),"\n"
