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
my $load_root = "/swnas/qc_nuclia";
my $proc_dir = "/project/PHG/PHG_Clinical/processing/$opt{prjid}";
unless ($opt{dir}) {
    $proc_dir = $opt{dir};
}
my $seqdatadir = "/project/PHG/PHG_Clinical/illumina";

$rinfo{'run.name'} = $opt{prjid};

@designfiles = `ls $proc_dir\/*/design.txt`;
chomp(@designfiles);

my $rnaseqnoumi;
if (-e "$seqdatadir\/$opt{prjid}/noumi/Stats/ConversionStats.xml") {
    $rnaseqnoumi = 1;
}

foreach my $dfile (@designfiles) {
    open IN, "<$dfile" or die $!;
    my $d = (split(/\//,$dfile))[-1];
    my $baitpool = (split(/\//,$dfile))[-2];
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
	$sinfo{$hash{SampleID}}{'dmux.conversion.stats'} = "$load_root/demultiplexing/$opt{prjid}/Stats/ConversionStats.xml";	
	if ($rnaseqnoumi && $baitpool =~ m/rnaseq/) {
	    $sinfo{$hash{SampleID}}{'dmux.conversion.stats'} = "$load_root/demultiplexing/$opt{prjid}/noumi/Stats/ConversionStats.xml";
	}
	$sinfo{$hash{SampleID}}{'bait.pool'} = $baitpool;
	$sinfo{$hash{SampleID}}{'project.name'}=$hash{FamilyID};
	$sinfo{$hash{SampleID}}{'sample.name'}=$hash{SampleID};
	my @seqstats = `ssh answerbe\@198.215.54.71 "find $load_root\/seqanalysis//$opt{prjid}/analysis -type f -name $hash{SampleID}.sequence.stats.txt"`;
	chomp(@seqstats);
	$sinfo{$hash{SampleID}}{'sample.alignment'} = $seqstats[0];
	my @exoncov = `ssh answerbe\@198.215.54.71 "find $load_root\/seqanalysis//$opt{prjid}/analysis -type f -name $hash{SampleID}_exoncoverage.txt"`;
	chomp(@exoncov);
	if ($exoncov[0]) {
	    $sinfo{$hash{SampleID}}{'sample.coverage.raw'} = $exoncov[0];
	    $uniqcov = $exoncov[0];
	    $uniqcov =~ s/exoncoverage/exoncoverageuniq/;
	    $sinfo{$hash{SampleID}}{'sample.coverage.uniq'} = $uniqcov;
	}
	my @tdfs = `ssh answerbe\@198.215.54.71 "find $load_root\/seqanalysis/$opt{prjid}/analysis -type f -name \"$hash{SampleID}*.tdf\""`;
	
	chomp(@tdfs);
	my @rawtdfs = grep(/raw/,@tdfs);
	my @uniqtdfs = grep(/uniq/,@tdfs); 
	$sinfo{$hash{SampleID}}{'tdf.raw'} = shift @rawtdfs;
	$sinfo{$hash{SampleID}}{'tdf.uniq'} = shift @uniqtdfs;
	my @somstats = `ssh answerbe\@198.215.54.71 "find $load_root\/seqanalysis//$opt{prjid}/analysis -type f -name $hash{FamilyID}.sequence.stats.txt"`;
	chomp(@somstats);
	$sinfo{$hash{SampleID}}{'somatic.seq.stats'}=$somstats[0];
	if ($hash{FamilyID} =~ m/ROS1/) {
	    $sinfo{$hash{SampleID}}{'somatic.translocation'} = "$load_root\/seqanalysis/$opt{'prjid'}/analysis/$hash{FamilyID}\/$hash{SampleID}\/$hash{SampleID}\.translocations.answer.txt";
	}
	if ($hash{FamilyID} =~ m/GM12878/) {
	    my @snsp = `ssh answerbe\@198.215.54.71 "find $load_root\/seqanalysis//$opt{prjid}/analysis/$hash{FamilyID} -type f -name $hash{SampleID}.snsp.txt"`;
	    chomp(@snsp);
	    $sinfo{$hash{SampleID}}{'giab.snsp'} = $snsp[0];
	}
	
    }
}

open SH, ">$proc_dir\/$opt{prjid}.curlcommand.sh" or die $!;
print SH "nucliatoken=\$1\n";
my @prop = ('dmux.conversion.stats','bait.pool','project.name','sample.alignment','sample.coverage.raw','sample.coverage.uniq',
	    'sample.name','somatic.seq.stats','tdf.raw','tdf.uniq','giab.snsp',
	    'somatic.translocation');

foreach $sid (sort {$a cmp $b} keys %sinfo) {
    open OUT, ">$proc_dir/$sid\.properties" or die $!;
    foreach $rprop (keys %rinfo) {
	print OUT join("=",$rprop,$rinfo{$rprop}),"\n";
    }
    foreach $prop (@prop) {
	$sinfo{$sid}{$prop} = '' unless ($sinfo{$sid}{$prop});
	print OUT join("=",$prop,$sinfo{$sid}{$prop}),"\n";
    }
    close OUT;
    print SH (qq{curl "http://nuclia.biohpc.swmed.edu:8080/NuCLIAVault/addPipelineResultsWithProp?token=\$nucliatoken&propFilePath=$load_root\/seqanalysis/$opt{'prjid'}/$sid\.properties"\n});
} 
close SH;
