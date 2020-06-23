#!/usr/bin/perl -w
#run_casava.pl

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my %opt = ();
my $results = GetOptions (\%opt,'help|h','input|i=s','output|o=s','umi|u=s',
			  'fqout|f=s','processdir|d=s','seqrunid|p=s',
			  'outnf|n=s','panelsdir|t=s');

if ($opt{umi}) {
  open RNASS, ">$opt{umi}" or die $!;
  $rnaseqnoumi = 1;
}
open SS, "<$opt{input}" or die $!;
open SSOUT, ">$opt{output}" or die $!;

my %spairs;
my %stype;
my %samples;

my $outdir = $opt{outdir};
my $outnf = $opt{outnf};
while (my $line = <SS>){
  chomp($line);
  $line =~ s/\r//g;
  $line =~ s/ //g;
  $line =~ s/,+$//g;
  if ($line =~ m/^\[Data\]/) {
    if ($opt{umi}) {
      print SSOUT join("\n","[Settings]","ReverseComplement,0","Read2UMILength,8","TrimUMI,1"),"\n";
    }
    if ($opt{umi}) {
      print RNASS $line,"\n";
    }
    print SSOUT $line,"\n";
    $header = <SS>;
    $header =~ s/\r//g;
    chomp($header);
    $header =~ s/Sample_*/Sample_/g;
    print SSOUT $header,"\n";
    if ($opt{umi}) {
      print RNASS $header,"\n";
    }    
    my @colnames = split(/,/,$header);
    while (my $line = <SS>) {
      chomp($line);
      $line =~ s/\r//g;
      $line =~ s/ //g;
      $line =~ s/,+$//g;
      my @row = split(/,/,$line);
      my %hash;
      foreach my $j (0..$#row) {
	$hash{$colnames[$j]} = $row[$j];
      }
      if ($hash{Sample_Name} =~ m/_Lib/) {
	  $hash{Sample_Name} =~ s/_Lib.*//;
      }
      $hash{Sample_Name} =~ s/T_RNA_panelrnaseq-\d+-\d+/T_RNA_panelrnaseq/;
      $hash{Sample_Project} = $hash{Project} if $hash{Project};
      $hash{Sample_Project} =~ s/\s*$//g;
      $hash{Assay} = lc($hash{Assay});
      if ($hash{Sample_Name} =~ m/panel1385v2|rnaseq|pancancer|heme/) {
	  $hash{Assay} = 'panel1385v2' if ($hash{Sample_Name} =~ m/panel1385v2/);
	  $hash{Assay} = 'tspcrnaseq' if ($hash{Sample_Name} =~ m/panelrnaseq/);
	  $hash{Assay} = 'idtrnaseq' if ($hash{Sample_Name} =~ m/panelrnaseq\d+/i);
	  $hash{Assay} = 'wholernaseq' if ($hash{Sample_Name} =~ m/wholernaseq/);
	  $hash{Assay} = 'pancancer' if ($hash{Sample_Name} =~ m/pancancer\d*/i);
	  $hash{Assay} = 'heme' if ($hash{Sample_Name} =~ m/heme\d*/i);
      }else {
	  $hash{Assay} = 'panel1385v2' if ($hash{Assay} =~ m/panel1385v2/);
	  $hash{Assay} = 'tspcrnaseq' if ($hash{Assay} =~ m/panelrnaseq/);
	  $hash{Assay} = 'idtrnaseq' if ($hash{Assay} =~ m/panelrnaseq\d+/i);
	  $hash{Assay} = 'wholernaseq' if ($hash{Assay} =~ m/wholernaseq/);
	  $hash{Assay} = 'pancancer' if ($hash{Assay} =~ m/pancancer\d*/i);
	  $hash{Assay} = 'heme' if ($hash{Assay} =~ m/heme\d*/i);
      }
      unless ($hash{Class}) {
	$hash{Class} = 'tumor';
	$hash{Class} = 'normal' if ($hash{Sample_Name} =~ m/_N_/);
      }
      $hash{SubjectID} = $hash{Sample_Project};
      $hash{Sample_ID} = $hash{Sample_Name};
      $spairs{$hash{Assay}}{$hash{SubjectID}}{lc($hash{Class})}{$hash{Sample_Name}} = 1;
      $stype{$hash{Assay}}{$hash{Sample_Name}} = $hash{SubjectID};
      push @{$samples{$hash{SubjectID}}}, $hash{Sample_Name};
      my @newline;
      foreach my $j (0..$#row) {
	push @newline, $hash{$colnames[$j]};
      }
      if ($opt{umi} && $hash{Index_Name} !~ m/UMI/) {
	print RNASS join(",",@newline),"\n";
      }else {
	print SSOUT join(",",@newline),"\n";
      }
    }
  } else {
    if ($opt{umi}) {
      print RNASS $line,"\n";
    }
    print SSOUT $line,"\n";
  }
}
close SSOUT;

my %panelbed=("panel1385"=>"UTSW_V2_pancancer", "panel1385v2"=>"UTSW_V3_pancancer", "medexomeplus"=>"UTSW_medexomeplus", "heme"=>"UTSW_V4_heme", "pancancer"=>"UTSW_V4_pancancer", "idtrnaseq"=>"UTSW_V4_rnaseq", "solid"=>"UTSW_V4_solid");

open DES, ">$opt{processdir}\/subjects.txt" or die $!;
open CAS, ">$opt{processdir}\/lnfq.sh" or die $!;
print CAS "#!/bin/bash\n";
foreach $subjid (keys %samples)  {
    system(qq{mkdir -p $opt{processdir}/analysis/$subjid});
    my $inseqdir =  "$opt{fqout}/$opt{seqrunid}/$subjid";
    my $outseqdir = "$opt{processdir}/analysis/$subjid/fastq";
    print DES $subjid,"\n";
    foreach $samp (@{$samples{$subjid}}) {
	print CAS "ln -fs $inseqdir/$samp*_R1_*.fastq.gz $opt{processdir}/fastq/$samp\.R1.fastq.gz\n";
	print CAS "ln -fs $inseqdir/$samp*_R2_*.fastq.gz $opt{processdir}/fastq/$samp\.R2.fastq.gz\n";
    }
}
close DES;

my %inpairs;
foreach $dtype (keys %spairs ){
  system(qq{mkdir -p $opt{processdir}/$dtype});
  my %vars;
  if ($dtype !~ m/rna/) {
    $vars{capturedir} = "$opt{panelsdir}/$panelbed{$dtype}";
    $vars{capturebed} = "$vars{capturedir}/targetpanel.bed";
    if (-e "$vars{capturedir}/mutect2.pon.vcf.gz") {
      $vars{mutectpon} = "$vars{capturedir}/mutect2.pon.vcf.gz";
    }
  }
  $vars{output}="$opt{processdir}/analysis";
  $vars{input}="$opt{processdir}/fastq";
  foreach my $subjid (keys %{$spairs{$dtype}}) {
    my $pct = 0;
    foreach $tumorid (keys %{$spairs{$dtype}{$subjid}{tumor}}) {
      $caseid=$subjid;
      if ($pct > 0) {
	$caseid .= ".$pct";
      }
      $pct ++;
      my $normalid;
      if ($spairs{$dtype}{$subjid}{normal}) {
	@nids = keys %{$spairs{$dtype}{$subjid}{normal}};
	foreach $normalid (@nids) {
	  $normalid = shift @nids;
	  $inpairs{$normalid} = [$caseid,$tumorid,$normalid];
	  $inpairs{$tumorid} = [$caseid,$tumorid,$normalid];
	}
      }
      if ($subjid =~ m/GM12878/) {
	$vars{vcf} = join("_",$caseid,$opt{seqrunid}).'.dna.vcf.gz';
	$vars{sampid} = $tumorid;
	$vars{bamfile} = $tumorid."/".$tumorid.".bam";
      }
    }
  }
  open VAR, ">$opt{processdir}/$dtype/vars.sh" or die $!;
  print VAR "#!/bin/bash\n";
  foreach $key (keys %vars) {
    print VAR join("=",$key,$vars{$key}),"\n";
  }
  close VAR;
}

my $load_root = "/swnas/qc_nuclia";
my $proc_dir = "$opt{processdir}/$opt{seqrunid}";
my $seqrunid = $opt{seqrunid};
if ($opt{dir}) {
  $proc_dir = $opt{dir};
}
my @prop = ('run.name','dmux.conversion.stats','bait.pool','project.name','sample.alignment',
	    'sample.coverage.raw','sample.name','somatic.seq.stats','giab.snsp',
	    'somatic.translocation');

my $seqdatadir = "/project/PHG/PHG_Clinical/illumina";
foreach $dtype (keys %stype) {
  open SSOUT, ">$opt{processdir}\/$dtype\/design.txt" or die $!;
  print SSOUT join("\t","SampleID",'CaseID','TumorID','NormalID','FqR1','FqR2'),"\n";
  foreach $sampleid (sort {$a cmp $b} keys %{$stype{$dtype}}) {
    my %info;
    my $tumorid = $sampleid;
    my $normalid = '';
    $caseid = $stype{$dtype}{$sampleid};
    if ($inpairs{$sampleid}) {
      ($caseid,$tumorid,$normalid) = @{$inpairs{$sampleid}};
    }
    print SSOUT join("\t",$sampleid,$caseid,$tumorid,$normalid,
		     "$sampleid.R1.fastq.gz","$sampleid.R2.fastq.gz"),"\n";
    my $casedir="$load_root/seqanalysis/$seqrunid/analysis/$caseid";
    $info{'run.name'} = $seqrunid;
    $info{'dmux.conversion.stats'} = "$load_root/demultiplexing/$seqrunid/Stats/ConversionStats.xml";
    if ($rnaseqnoumi && $dtype =~ m/rnaseq/) {
      $info{'dmux.conversion.stats'} = "$load_root/demultiplexing/$seqrunid/noumi/Stats/ConversionStats.xml";
    }
    $info{'bait.pool'} = $dtype;
    $info{'project.name'}=$caseid;
    $info{'sample.name'}=$sampleid;
    $info{'sample.alignment'} = "$casedir/$sampleid/$sampleid.sequence.stats.txt";
    $info{'sample.coverage.raw'} = "$casedir/$sampleid/$sampleid\_exoncoverage.txt";
    if ($inpairs{$sampleid}) {
      $info{'somatic.seq.stats'}="$casedir/dna_$seqrunid/$caseid\_$seqrunid.sequence.stats.txt";
    }
    if ($caseid =~ m/ROS1/) {
      $info{'somatic.translocation'} = "$casedir/$sampleid/$sampleid.translocations.answer.txt";
    }
    if ($caseid =~ m/GM12878/) {
      $info{'giab.snsp'} = "$casedir/$sampleid.snsp.txt";
    }
    open OUT, ">$opt{processdir}/$sampleid.properties" or die $!;
    foreach $prop (@prop) {
      $info{$prop} = '' unless ($info{$prop});
      print OUT join("=",$prop,$info{$prop}),"\n";
    }
  }
}
