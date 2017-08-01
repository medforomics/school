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
      $clinres = 'toresearch' if ($hash{Description} && $hash{Description} =~ m/research/i);
      $hash{ClinRes} = $clinres;
      $hash{Sample_Project} = $hash{Project} if $hash{Project};
      $hash{Sample_Project} =~ s/\s*$//g;
      my @samplename = split(/_/,$hash{Sample_Name});
      unless ($hash{Class}) {
	$hash{Class} = 'tumor';
	$hash{Class} = 'normal' if ($hash{Sample_Name} =~ m/_N_/);
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
      $samps{$hash{SubjectID}}{lc($hash{Class})} = $hash{MergeName} if ($hash{Assay} =~ m/panel1385|medexome/i);
      $sampleinfo{$samp} = \%hash;
      push @{$samples{lc($hash{Assay})}{$hash{SubjectID}}}, $samp;
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
print CAS "#SBATCH -t 14-0:0:00\n#SBATCH -o $prjid.out\n#SBATCH -e $prjid.err\n";
print CAS "#SBATCH --mail-type ALL\n#SBATCH --mail-user erika.villa\@utsouthwestern.edu\n";
print CAS "source /etc/profile.d/modules.sh\n";
print CAS "module load bcl2fastq/2.17.1.14 fastqc/0.11.2 nextflow/0.24.1-SNAPSHOT\n";
my $seqdatadir = "/project/PHG/PHG_Illumina/BioCenter/$prjid";
if (-e "/project/PHG/PHG_Illumina/Research/$prjid") {
  $seqdatadir = "/project/PHG/PHG_Illumina/Research/$prjid";
}

print CAS "bcl2fastq --barcode-mismatches 0 -o /project/PHG/PHG_Clinical/illumina/$prjid --ignore-missing-positions --no-lane-splitting --ignore-missing-filter --ignore-missing-bcls --runfolder-dir $seqdatadir --sample-sheet /project/PHG/PHG_Clinical/illumina/sample_sheets/$prjid\.bcl2fastq.csv &> /project/PHG/PHG_Clinical/illumina/logs/run_casava_$prjid\.log\n";
print CAS "mkdir /project/PHG/PHG_BarTender/bioinformatics/demultiplexing/$prjid\n" unless (-e "/project/PHG/PHG_BarTender/bioinformatics/demultiplexing/$prjid");
print CAS "mv /project/PHG/PHG_Clinical/illumina/$prjid\/Reports /project/PHG/PHG_BarTender/bioinformatics/demultiplexing/$prjid\n" unless (-e "/project/PHG/PHG_BarTender/bioinformatics/demultiplexing/$prjid/Reports");
print CAS "mv /project/PHG/PHG_Clinical/illumina/$prjid\/Stats /project/PHG/PHG_BarTender/bioinformatics/demultiplexing/$prjid\n" unless (-e "/project/PHG/PHG_BarTender/bioinformatics/demultiplexing/$prjid/Stats");

foreach $dtype (keys %samples) {
  my $prodir = "/project/PHG/PHG_Clinical/processing";
  my $outdir = "$prodir\/$prjid/fastq";
  my $outnf = "$prodir\/$prjid/analysis";
  my $workdir = "$prodir\/$prjid/work";
  system("mkdir $prodir\/$prjid") unless (-e "$prodir\/$prjid");
  system("mkdir $outdir") unless (-e $outdir);
  system("mkdir $outnf") unless (-e $outnf);
  system("mkdir $workdir") unless (-e $workdir);
  my %completeout; 
  open SSOUT, ">$outdir\/design.txt" or die $!;
  print SSOUT join("\t","SampleMergeName",'SampleID','SampleName','SubjectID','FullPathToFqR1','FullPathToFqR2'),"\n";
  my @copyfq;
  foreach $project (keys %{$samples{$dtype}}) {
    my $datadir =  "/project/PHG/PHG_Clinical/illumina/$prjid/$project/";
    foreach $samp (@{$samples{$dtype}{$project}}) {
      my %info = %{$sampleinfo{$samp}};
      print CAS "ln -s $datadir/$samp*_R1_*.fastq.gz $outdir\/$samp\.R1.fastq.gz\n";
      print CAS "ln -s $datadir/$samp*_R2_*.fastq.gz $outdir\/$samp\.R2.fastq.gz\n";
      my $finaloutput = '/project/PHG/PHG_Clinical/'.$info{ClinRes};
      unless (-e "$finaloutput\/$info{SubjectID}") {
	  system("mkdir $finaloutput\/$info{SubjectID}");
      }
      my $finalrestingplace = "$finaloutput\/$info{SubjectID}\/$info{MergeName}";
      if (-e $finalrestingplace) {
	  $finalrestingplace .= "_".(split(/_|-/,$prjid))[-1];
      }
      system("mkdir $finalrestingplace");
      $completeout{$info{MergeName}} = $finalrestingplace;
      print CAS "ln -s $datadir/$samp*_R1_*.fastq.gz $finalrestingplace\/$samp\.R1.fastq.gz\n";
      print CAS "ln -s $datadir/$samp*_R2_*.fastq.gz $finalrestingplace\/$samp\.R2.fastq.gz\n";
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
      foreach $clast (@ctypes) {
	print TONLY join("\t",$samps{$subjid}{$clast},
			 $samps{$subjid}{$clast}.".bam",
			 $samps{$subjid}{$clast}.".final.bam"),"\n";
	$tonlys ++;
      }
      if ($samps{$subjid}{tumor} && $samps{$subjid}{normal}) {
	print TNPAIR join("\t",$samps{$subjid}{tumor},$samps{$subjid}{normal},
			  $samps{$subjid}{tumor}.".bam",
			  $samps{$subjid}{normal}.".bam",
			  $samps{$subjid}{tumor}.".final.bam",
			  $samps{$subjid}{normal}.".final.bam"),"\n";
	$tnpairs ++;
      }
    }
    close TNPAIR;
    close TONLY;
  }
  my $capture = '/project/shared/bicf_workflow_ref/GRCh38/UTSWV2.bed';
  $capture = '/project/shared/bicf_workflow_ref/GRCh38/MedExome_Plus.bed' if ($dtype eq 'medexomeplus');
  print CAS "cd /project/PHG/PHG_Clinical/processing/$prjid\n";
  if ($dtype eq 'panel1385' || $dtype eq 'medexomeplus') {
    print CAS "nextflow -C /project/PHG/PHG_Clinical/clinseq_workflows/nextflow.config run -w $workdir /project/PHG/PHG_Clinical/clinseq_workflows/alignment.nf --design $outdir\/design.txt --capture $capture --input $outdir --output $outnf >> $outnf\/nextflow_alignment.log\n";
    print CAS "nextflow -C /project/PHG/PHG_Clinical/clinseq_workflows/nextflow.config run -w $workdir /project/PHG/PHG_Clinical/clinseq_workflows/somatic.nf --design $outdir\/design_tumor_normal.txt --capture $capture --input $outnf --output $outnf >> $outnf\/nextflow_somatic.log &\n" if ($tnpairs);
    print CAS "nextflow -C /project/PHG/PHG_Clinical/clinseq_workflows/nextflow.config run -w $workdir /project/PHG/PHG_Clinical/clinseq_workflows/tumoronly.nf --design $outdir\/design_tumor_only.txt --capture $capture --input $outnf --output $outnf >> $outnf\/nextflow_tumoronly.log &\n" if ($tonlys);
    print CAS "wait\n";
  }elsif ($dtype eq 'panelrnaseq' || $dtype eq 'wholernaseq') {
    my $mdup = 'skip';
    $mdup = 'mark' if ($dtype eq 'wholernaseq');
    print CAS "nextflow -C /project/PHG/PHG_Clinical/clinseq_workflows/nextflow.config run -w $workdir /project/PHG/PHG_Clinical/clinseq_workflows/rnaseq.nf --design $outdir\/design.txt --input $outdir --output $outnf --markdups $mdup &> $outnf\/nextflow_rnaseq.log &\n";
    print CAS "wait\n";
  }
  foreach $sampid (keys %completeout) {
      $finalrestingplace = $completeout{$sampid};
      print CAS "mv $outnf\/$sampid\* $finalrestingplace\n";
  }
}
close CAS;
