#!/usr/bin/perl -w
#create_properties_run.pl

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my %opt = ();
my $results = GetOptions (\%opt,'help|h','prjid|p=s','dir|d=s');

if (!defined $opt{prjid} || $opt{help}) {
  $usage = <<EOF;
  usage: $0 -p prjid
      
      -p prjid -- this is the project name in /project/PHG/PHG_Illumina/BioCenter/ 140505_SN7001189_0117_AH7LRLADXX
      
EOF
  die $usage,"\n";
}

my %rinfo;
my %sinfo;

$rinfo{'dmux.conversion.stats'} = "/project/PHG/PHG_BarTender/bioinformatics/demultiplexing/$opt{prjid}/Stats/ConversionStats.xml";
$rinfo{'run.name'} = $opt{prjid};

my @designfiles = `ls $opt{dir}/$opt{prjid}/fastq/*.design.txt`;
unless (@designfiles) {
    @designfiles = `ls $opt{dir}/$opt{prjid}/*/design.txt`;
}

chomp(@designfiles);
foreach my $dfile (@designfiles) {
    open IN, "<$dfile" or die $!;
    my $d = (split(/\//,$dfile))[-1];
    my $baitpool;
    if ($d eq 'design.txt') {
	$baitpool = (split(/\//,$dfile))[-2];
    }else {
	$baitpool = (split(/\./,$d))[0];
    }
    next if ($baitpool eq 'wholernaseq');
    my $header = <IN>;
    chomp($header);
    my @colnames = split(/\t/,$header);
    while (my $line = <IN>) {
	chomp($line);
	my @row = split(/\t/,$line);
	my %hash;
	foreach my $i (0..$#row) {
	    $hash{$colnames[$i]} = $row[$i];
	}
	$sinfo{$hash{SampleID}}{'bait.pool'} = $baitpool;
	$sinfo{$hash{SampleID}}{'project.name'}=$hash{FamilyID};
	$sinfo{$hash{SampleID}}{'sample.name'}=$hash{SampleID};
	my @seqstats = `find $opt{dir}/$opt{prjid}/analysis -type f -name $hash{SampleID}.sequence.stats.txt`;
	chomp(@seqstats);
	$sinfo{$hash{SampleID}}{'sample.alignment'} = $seqstats[0];
	my @exoncov = `find $opt{dir}/$opt{prjid}/analysis -name $hash{SampleID}_exoncoverage.txt`;
	chomp(@exoncov);
	if ($exoncov[0]) {
	    $sinfo{$hash{SampleID}}{'sample.coverage.raw'} = $exoncov[0];
	    $uniqcov = $exoncov[0];
	    $uniqcov =~ s/exoncoverage/exoncoverageuniq/;
	    $sinfo{$hash{SampleID}}{'sample.coverage.uniq'} = $uniqcov;
	}
	my @tdfs = `find $opt{dir}/$opt{prjid}/analysis -type f -name $hash{SampleID}*.tdf`;
	chomp(@tdfs);
	@rtdfs = grep(/raw/,@tdfs);
	@utdfs = grep(/uniq/,@tdfs); 
	$sinfo{$hash{SampleID}}{'tdf.raw'} = $rtdfs[0];
	$sinfo{$hash{SampleID}}{'tdf.uniq'} = $utdfs[0];
	my @somstats = `find $opt{dir}/$opt{prjid}/analysis -type f -name $hash{FamilyID}*sequence.stats.txt`;
	chomp(@somstats);
	$sinfo{$hash{SampleID}}{'somatic.seq.stats'}=$somstats[0];
	if ($hash{FamilyID} =~ m/ROS1/) {
	    $sinfo{$hash{SampleID}}{'somatic.translocation'} = "$opt{dir}\/$opt{'prjid'}/analysis/$hash{FamilyID}\/$hash{SampleID}\/$hash{SampleID}\.translocations.txt";
	}
	if ($hash{FamilyID} =~ m/GM12878/) {
	    my @snsp = `find $opt{dir}\/$opt{prjid}\/ -type f -name "*.snsp.txt"`;
	    chomp(@snsp);
	    $sinfo{$hash{SampleID}}{'giab.snsp'} = $snsp[0];
	}
	
    }
}
open SH, ">$opt{dir}\/$opt{prjid}.curlcommand.sh" or die $!;
print SH "nucliatoken=\$1\n";
my @prop = ('bait.pool','project.name','sample.alignment','sample.coverage.raw','sample.coverage.uniq',
	    'sample.name','somatic.seq.stats','tdf.raw','tdf.uniq','giab.snsp',
	    'somatic.translocation');

foreach $sid (sort {$a cmp $b} keys %sinfo) {
    open OUT, ">$opt{dir}\/$opt{prjid}/$sid\.properties" or die $!;
    print OUT join("=",'dmux.conversion.stats',$rinfo{'dmux.conversion.stats'}),"\n";
    print OUT join("=",'run.name',$rinfo{'run.name'}),"\n";
    foreach $prop (@prop) {
	$sinfo{$sid}{$prop} = '' unless ($sinfo{$sid}{$prop});
	print OUT join("=",$prop,$sinfo{$sid}{$prop}),"\n";
    }
    close OUT;
    $outdir = $opt{dir};
    print SH (qq{curl "http://nuclia.biohpc.swmed.edu:8080/NuCLIAVault/addPipelineResultsWithProp?token=\$nucliatoken&propFilePath=$outdir\/$opt{prjid}/$sid\.properties"\n});
} 
close SH;
