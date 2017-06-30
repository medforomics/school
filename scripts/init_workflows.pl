#!/usr/bin/perl -w
#run_casava.pl

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my %opt = ();
my $results = GetOptions (\%opt,'help|h','prjid|p=s');

if (!defined $opt{prjid} || $opt{help}) {
  $usage = <<EOF;
usage: run_casava.pl -p prjid

-p prjid -- this is the project name in /project/PHG/PHG_Illumina/BioCenter/ 140505_SN7001189_0117_AH7LRLADXX

EOF
  die $usage,"\n";
}

my $prjid = $opt{prjid};

open SS, "</project/PHG/PHG_Clinical/illumina/sample_sheets/$prjid\.csv" or die $!;
open SSOUT, ">/project/PHG/PHG_Clinical/illumina/sample_sheets/$prjid\.bcl2fastq.csv" or die $!;
my %sampleinfo;
while (my $line = <SS>){
  chomp($line);
  $line =~ s/\r//g;
  $line =~ s/ //g;
  $line =~ s/,+$//g;
  print SSOUT $line,"\n";
  if ($line =~ m/^\[Data\]/) {
    $header = <SS>;
    $header =~ s/\r//g;
    chomp($header);
    $header =~ s/Sample_*/Sample_/g;
    print SSOUT $header,"\n";
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
      my $clinres = 'complete';
      $clinres = 'toresearch' if ($hash{Description} =~ m/research/i);
      $hash{Sample_Project} = $hash{Project} if $hash{Project};
      $hash{Sample_Project} =~ s/\s*$//g;
      $project = $hash{Sample_Project};
      my @samplename = split(/_/,$hash{Sample_Name});
      unless ($hash{Class}) {
	$hash{Class} = 'Tumor';
	$hash{Class} = 'Normal' if ($hash{Sample_Name} =~ m/_N_/);
      }
      $hash{SubjectID} = $hash{Sample_Project};
      unless ($hash{MergeName}) {
	$hash{MergeName} = $hash{Sample_Name};
	if ($samplename[-1] =~ m/^[A|B|C|D]$/) {
	  pop @samplename;
	  $hash{MergeName} = join("_",@samplename);
	}
      }
      $hash{Sample_Name} = join("_ClarityID-",$hash{Sample_Name},$hash{Sample_ID});
      $hash{Sample_ID} = $hash{Sample_Name};
      $samp = $hash{Sample_Name};
      $dtype = lc($hash{Assay});
      $subjgrp{$hash{SubjectID}}{$hash{MergeName}} = $clinres;
      $samps{$hash{SubjectID}}{lc($hash{Class})} = $hash{MergeName} if ($dtype =~ m/panel1385|medexome/);
      $sampleinfo{$samp} = \%hash;
      push @{$samples{$dtype}{$project}}, $samp;
      my @newline;
      foreach my $j (0..$#row) {
	push @newline, $hash{$colnames[$j]};
      }
      print SSOUT join(",",@newline),"\n";
    }
  }
}
close SSOUT;

open CAS, ">/project/PHG/PHG_Clinical/illumina/logs/run_casava_$prjid\.sh" or die $!;
print CAS "#!/bin/bash\n#SBATCH --job-name $prjid\n#SBATCH -N 1\n";
print CAS "#SBATCH -t 3-0:0:00\n#SBATCH -o $prjid.out\n#SBATCH -e $prjid.err\n";
print CAS "#SBATCH --mail-type ALL\n#SBATCH --mail-user erika.villa\@utsouthwestern.edu\n";
print CAS "source /etc/profile.d/modules.sh\n";
print CAS "module load bcl2fastq/2.17.1.14 fastqc/0.11.2 nextflow/0.24.1-SNAPSHOT\n";
print CAS 'export PERL5LIB=/project/BICF/BICF_Core/shared/seqprg/lib/perl5/x86_64-linux-thread-multi/:$PERL5LIB',"\n";
my $seqdatadir = "/project/PHG/PHG_Illumina/BioCenter/$prjid";
if (-e "/project/PHG/PHG_Illumina/Research/$prjid") {
  $seqdatadir = "/project/PHG/PHG_Illumina/Research/$prjid";
}

print CAS "bcl2fastq --barcode-mismatches 0 -o /project/PHG/PHG_Clinical/illumina/$prjid --ignore-missing-positions --no-lane-splitting --ignore-missing-filter --ignore-missing-bcls --runfolder-dir $seqdatadir --sample-sheet /project/PHG/PHG_Clinical/illumina/sample_sheets/$prjid\.bcl2fastq.csv &> /project/PHG/PHG_Clinical/illumina/logs/run_casava_$prjid\.log\n";
print CAS "mkdir /project/PHG/PHG_BarTender/bioinformatics/demultiplexing/$prjid\n" unless (-e "/project/PHG/PHG_BarTender/bioinformatics/demultiplexing/$prjid");
print CAS "mv /project/PHG/PHG_Clinical/illumina/$prjid\/Reports /project/PHG/PHG_BarTender/bioinformatics/demultiplexing/$prjid\n" unless (-e "/project/PHG/PHG_BarTender/bioinformatics/demultiplexing/$prjid/Reports");
print CAS "mv /project/PHG/PHG_Clinical/illumina/$prjid\/Stats /project/PHG/PHG_BarTender/bioinformatics/demultiplexing/$prjid\n" unless (-e "/project/PHG/PHG_BarTender/bioinformatics/demultiplexing/$prjid/Stats");

foreach $dtype (keys %samples) {
  my $outdir = "/project/PHG/PHG_Clinical/processing/$prjid/fastq";
  my $outnf = "/project/PHG/PHG_Clinical/processing/$prjid/analysis";
  my $workdir = "/project/PHG/PHG_Clinical/processing/$prjid/work";
  system("mkdir /project/PHG/PHG_Clinical/processing/$prjid");
  system("mkdir $outdir");
  system("mkdir $outnf");
  system("mkdir $workdir");
  open SSOUT, ">$outdir\/design.txt" or die $!;
  print SSOUT join("\t","SampleMergeName",'SampleID','SampleName','SubjectID','FullPathToFqR1','FullPathToFqR2'),"\n";
  foreach $project (keys %{$samples{$dtype}}) {
    my $datadir =  "/project/PHG/PHG_Clinical/illumina/$prjid/$project/";
    foreach $samp (@{$samples{$dtype}{$project}}) {
      my %info = %{$sampleinfo{$samp}};
      print CAS "ln -s $datadir/$samp*_R1_*.fastq.gz $outdir\/$samp\.R1.fastq.gz\n";
      print CAS "ln -s $datadir/$samp*_R2_*.fastq.gz $outdir\/$samp\.R2.fastq.gz\n";
      print SSOUT join("\t",$info{MergeName},$info{Sample_ID},$info{Sample_Name},$info{SubjectID},"$samp\.R1.fastq.gz","$samp\.R2.fastq.gz"),"\n";
    }
  }
  close SSOUT;
  my $tnpairs = 0;
  my $tonlys = 0;
  if ($dtype =~ m/panel1385|medexome/) {
    open TONLY, ">$outdir\/design_tumor_only.txt" or die $!;
    open TNPAIR, ">$outdir\/design_tumor_normal.txt" or die $!;
    print TNPAIR join("\t",'TumorID','NormalID','TumorBAM','NormalBAM',
		      'TumorOntargetBAM','NormalOntargetBAM'),"\n";
    print TONLY join("\t",'SampleID','BAM','OntargetBAM'),"\n";
    foreach my $subjid (keys %samps) {
      my @ctypes = keys %{$samps{$subjid}};
      if (scalar(@ctypes) > 1) {
	print TNPAIR join("\t",$samps{$subjid}{tumor},$samps{$subjid}{normal},
			  $samps{$subjid}{tumor}.".bam",$samps{$subjid}{normal}.".bam",
			  $samps{$subjid}{tumor}.".final.bam",
			  $samps{$subjid}{normal}.".final.bam"),"\n";
	$tnpairs ++;
      }
      print TONLY join("\t",$samps{$subjid}{$ctypes[0]},$samps{$subjid}{$ctypes[0]}.".bam",
		       $samps{$subjid}{$ctypes[0]}.".final.bam"),"\n";
      $tonlys ++;
    }
    close TNPAIR;
    close TONLY;
  }
  my $capture = '/project/shared/bicf_workflow_ref/GRCh38/UTSWV2.bed';
  $capture = '/project/shared/bicf_workflow_ref/GRCh38/MedExome_Plus.bed' if ($dtype eq 'medexomeplus');
  print CAS "cd /project/PHG/PHG_Clinical/processing/$prjid\n";
  if ($dtype eq 'panel1385' || $dtype eq 'medexomeplus') {
    print CAS "nextflow -C /project/PHG/PHG_Clinical/clinseq_workflows/nextflow.config run -with-timeline -w $workdir /project/PHG/PHG_Clinical/clinseq_workflows/alignment.nf --design $outdir\/design.txt --capture $capture --input $outdir --output $outnf >> $outnf\/nextflow_alignment.log\n";
    print CAS "nextflow -C /project/PHG/PHG_Clinical/clinseq_workflows/nextflow.config run -with-timeline -w $workdir /project/PHG/PHG_Clinical/clinseq_workflows/somatic.nf --design $outdir\/design_tumor_normal.txt --capture $capture --input $outnf --output $outnf >> $outnf\/nextflow_somatic.log &\n" if ($tnpairs);
    print CAS "nextflow -C /project/PHG/PHG_Clinical/clinseq_workflows/nextflow.config run -with-timeline -w $workdir /project/PHG/PHG_Clinical/clinseq_workflows/tumoronly.nf --design $outdir\/design_tumor_only.txt --capture $capture --input $outnf --output $outnf >> $outnf\/nextflow_tumoronly.log &\n" if ($tonlys);
    print CAS "wait\n";
  }elsif ($dtype eq 'panelrnaseq' || $dtype eq 'wholernaseq') {
    my $mdup = 'skip';
    $mdup = 'mark' if ($dtype eq 'wholernaseq');
    print CAS "nextflow -C /project/PHG/PHG_Clinical/clinseq_workflows/nextflow.config run -with-timeline -w $workdir /project/PHG/PHG_Clinical/clinseq_workflows/rnaseq.nf --design $outdir\/design.txt --input $outdir --output $outnf --markdups $mdup &> $outnf\/nextflow_rnaseq.log &\n";
    print CAS "wait\n";
  }
  foreach $subjacc (keys %subjgrp) {
    my @samples = keys %{$subjgrp{$subjacc}};
    foreach $sampid (@samples) {
      my $finaloutput = '/project/PHG/PHG_Clinical/'.$subjgrp{$subjacc}{$sampid};
      unless (-e "$finaloutput\/$subjacc") {
	system("mkdir $finaloutput\/$subjacc");
      }
      my $finalrestingplace = "$finaloutput\/$subjacc\/$sampid";
      if (-e "$finaloutput\/$subjacc\/$sampid") {
	$finalrestingplace = "$finaloutput\/$subjacc\/$sampid".(split(/_/,$prjid))[0];
      }
      print CAS "mkdir $finalrestingplace\n";
      print CAS "mv $outnf\/$sampid\* $finalrestingplace\n";
    }
  }
}
close CAS;
