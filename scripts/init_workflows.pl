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
my $finaloutput = '/project/PHG/PHG_Clinical/validation';


open SS, "</project/PHG/PHG_Illumina/sample_sheets/$prjid\.csv" or die $!;
open SSOUT, ">/project/PHG/PHG_Illumina/sample_sheets/$prjid\.bcl2fastq.csv" or die $!;
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
      $hash{Sample_Project} = $hash{Project} if $hash{Project};
      $hash{Sample_Project} =~ s/\s*$//g;
      $project = $hash{Sample_Project};
      $hash{Sample_Name} =~ s/_1385_/_panel1385_/;
      $hash{Sample_Name} =~ s/_N_/-N_/;
      $hash{Sample_Name} =~ s/_T_/-T_/;
      $hash{Sample_Name} =~ s/Hyb_/Hyb/;
      $hash{Sample_Name} =~ s/^1793/S1793/;
      my @samplename = split(/_/,$hash{Sample_Name});
      unless ($hash{Class}) {
	  $hash{Class} = 'Tumor';
	  $hash{Class} = 'Normal' if ($hash{Sample_Name} =~ m/-N_/);
      }if ($hash{Sample_Name} =~ m/$hash{Sample_Project}/) {
	  $hash{SubjectID} = $hash{Sample_Project};
      }else {
	  $hash{SubjectID} = $samplename[0];
	  $hash{SubjectID} =~ s/-N//;
	  $hash{SubjectID} =~ s/-T//;      
      }
      unless ($hash{SampleMergeName}) {
	  $hash{SampleMergeName} = $hash{Sample_Name};
	  if ($samplename[-1] =~ m/^[A|B|C|D]$/) {
	      pop @samplename;
	      $hash{SampleMergeName} = join("_",@samplename);
	  }
      }
      $hash{Sample_ID} = $hash{Sample_Name};
      $samp = $hash{Sample_Name};	
      $dtype = lc($hash{Assay});
      
      push @{$subjgrp{$hash{SubjectID}}}, $hash{SampleMergeName};
      $samps{$hash{SubjectID}}{lc($hash{Class})} = $hash{SampleMergeName} if ($dtype =~ m/panel1385|medexome/);
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

open CAS, ">/project/PHG/PHG_Illumina/logs/run_casava_$prjid\.sh" or die $!;
print CAS "#!/bin/bash\n#SBATCH --job-name bcl2fastq\n#SBATCH -N 1\n";
print CAS "#SBATCH -t 1-0:0:00\n#SBATCH -o job_%j.out\n#SBATCH -e job_%j.err\n";
print CAS "#SBATCH --mail-type ALL\n#SBATCH --mail-user erika.villa\@utsouthwestern.edu\n";
print CAS "module load bcl2fastq/2.17.1.14 fastqc/0.11.2 nextflow/0.20.1\n";
print CAS 'export PERL5LIB=/project/BICF/BICF_Core/shared/seqprg/lib/perl5/x86_64-linux-thread-multi/:$PERL5LIB',"\n";
my $seqdatadir = "/project/PHG/PHG_Illumina/BioCenter/$prjid";
if (-e "/project/PHG/PHG_Illumina/Research/$prjid") {
    $seqdatadir = "/project/PHG/PHG_Illumina/Research/$prjid";
}

print CAS "bcl2fastq --barcode-mismatches 0 -o /project/PHG/PHG_Clinical/illumina/$prjid --ignore-missing-positions --no-lane-splitting --ignore-missing-filter --ignore-missing-bcls --runfolder-dir $seqdatadir --sample-sheet /project/PHG/PHG_Illumina/sample_sheets/$prjid\.bcl2fastq.csv &> /project/PHG/PHG_Illumina/logs/run_casava_$prjid\.log\n";
#print CAS "/project/PHG/PHG_Illumina/scripts/parse_conversion_stats.pl /project/PHG/PHG_Clinical/illumina/$prjid\n";
foreach $dtype (keys %samples) {
  my $outdir = "/project/PHG/PHG_Clinical\/$dtype\/$prjid/fastq";
  my $outnf = "/project/PHG/PHG_Clinical\/$dtype\/$prjid/analysis";
  my $workdir = "/project/PHG/PHG_Clinical\/$dtype\/$prjid/work";
  system("mkdir /project/PHG/PHG_Clinical\/$dtype\/$prjid");
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
      print SSOUT join("\t",$info{SampleMergeName},$info{Sample_ID},$info{Sample_Name},$info{SubjectID},
		       "$samp\.R1.fastq.gz","$samp\.R2.fastq.gz"),"\n";
    }
  }
  close SSOUT;
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
			  $samps{$subjid}{tumor}.".ontarget.bam",
			  $samps{$subjid}{normal}.".ontarget.bam"),"\n";
      } else {
	print TONLY join("\t",$samps{$subjid}{$ctypes[0]},$samps{$subjid}{$ctypes[0]}.".bam",
			 $samps{$subjid}{$ctypes[0]}.".ontarget.bam"),"\n";
	  }
    }
    close TNPAIR;
    close TONLY;
  }
  my $capture = '/project/shared/bicf_workflow_ref/GRCh38/UTSWV2.bed';
  $capture = '/project/shared/bicf_workflow_ref/GRCh38/MedExome_Plus.bed' if ($dtype eq 'medexomeplus');
  print CAS "cd /project/PHG/PHG_Clinical\/$dtype\/$prjid\n";
  if ($dtype eq 'panel1385' || $dtype eq 'medexomeplus') {
    print CAS "nextflow -C /project/PHG/PHG_Clinical/clinseq_workflows/nextflow.config run -with-timeline -w $workdir /project/PHG/PHG_Clinical/clinseq_workflows/alignment.nf --design $outdir\/design.txt --capture $capture --input $outdir --output $outnf &> $outnf\/nextflow.log\n";
    print CAS "nextflow -C /project/PHG/PHG_Clinical/clinseq_workflows/nextflow.config run -with-timeline -w $workdir /project/PHG/PHG_Clinical/clinseq_workflows/somatic.nf --design $outdir\/design_tumor_normal.txt --capture $capture --input $outnf --output $outnf &> $outnf\/nextflow.log\n &";
    print CAS "nextflow -C /project/PHG/PHG_Clinical/clinseq_workflows/nextflow.config run -with-timeline -w $workdir /project/PHG/PHG_Clinical/clinseq_workflows/tumoronly.nf --design $outdir\/design_tumor_only.txt --capture $capture --input $outnf --output $outnf &> $outnf\/nextflow.log\n &";
    print CAS "wait\n";
  }elsif ($dtype eq 'panelrnaseq' || $dtype eq 'wholernaseq') {
    my $mdup = 'skip';
    $mdup = 'mark' if ($dtype eq 'wholernaseq');
    print CAS "nextflow -C /project/PHG/PHG_Clinical/clinseq_workflows/nextflow.config run -with-timeline -w $workdir /project/PHG/PHG_Clinical/clinseq_workflows/rnaseq.nf --design $outdir\/design.txt --input $outdir --output $outnf --markdups $mdup &> $outnf\/nextflow.log &\n";
    print CAS "wait\n";
  }
  foreach $subjacc (keys %subjgrp) {
      my @samples = @{$subjgrp{$subjacc}}
      unless (-e "$finaloutput\/$subjacc") {
	  print CAS "mkdir $finalrestingplace\n";
      }
      foreach $sampid (@samples) {
	  my $finalrestingplace = "$finaloutput\/$subjacc\/$sampid";
	  if (-e "$finaloutput\/$subjacc\/$sampid") {
	      $finalrestingplace = "$finaloutput\/$subjacc\/$sampid_".(split(/_/,$prjid))[0];
	  }
	  print CAS "mkdir $finalrestingplace\n";
	  print CAS "mv $outnf\/$sampid\.* $outnf\/$sampid\_* $finalrestingplace\n";
      }
  }
}
close CAS;
