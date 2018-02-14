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
my $oriss = "/project/PHG/PHG_Clinical/illumina/sample_sheets/$prjid\.csv";
my $newss = "/project/PHG/PHG_Clinical/illumina/sample_sheets/$prjid\.bcl2fastq.csv";

my $seqdatadir = "/project/PHG/PHG_Illumina/BioCenter/$prjid";
if (-e "/project/PHG/PHG_Illumina/Research/$prjid") {
  $seqdatadir = "/project/PHG/PHG_Illumina/Research/$prjid";
}

$umi = `grep "<Read Number=\\\"2\\\" NumCycles=\\\"14\\\" IsIndexedRead=\\\"Y\\\" />" $seqdatadir/RunInfo.xml`;

open SS, "<$oriss" or die $!;
open SSOUT, ">$newss" or die $!;

my %sampleinfo;
while (my $line = <SS>){
  chomp($line);
  $line =~ s/\r//g;
  $line =~ s/ //g;
  $line =~ s/,+$//g;
  if ($line =~ m/^\[Data\]/) {
    if ($umi) {
      print SSOUT join("\n","[Settings]","ReverseComplement,0","Read2UMILength,8"),"\n";
    }
    print SSOUT $line,"\n";
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
      $hash{Assay} = 'panel1385' if ($hash{Assay} eq 'dnaseqdevelopment');
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
      if ($samps{$hash{SubjectID}}{lc($hash{Class})} && 
	  $samps{$hash{SubjectID}}{lc($hash{Class})} ne $hash{MergeName}) {
	  $hash{Class} = join('.',$hash{Class},$.);
      }
      $samps{$hash{SubjectID}}{lc($hash{Class})} = $hash{MergeName};
      $sampleinfo{$samp} = \%hash;
      push @{$samples{lc($hash{Assay})}{$hash{SubjectID}}}, $samp;
      my @newline;
      foreach my $j (0..$#row) {
	  push @newline, $hash{$colnames[$j]};
      }
      print SSOUT join(",",@newline),"\n";
    }
  } else {
    print SSOUT $line,"\n";
  }
}
close SSOUT;

open CAS, ">/project/PHG/PHG_Clinical/illumina/logs/run_casava_$prjid\.sh" or die $!;
print CAS "#!/bin/bash\n#SBATCH --job-name $prjid\n#SBATCH -N 1\n";
print CAS "#SBATCH -t 14-0:0:00\n#SBATCH -o $prjid.out\n#SBATCH -e $prjid.err\n";
print CAS "#SBATCH --mail-type ALL\n#SBATCH --mail-user erika.villa\@utsouthwestern.edu\n";
print CAS "source /etc/profile.d/modules.sh\n";
print CAS "module load bcl2fastq/2.17.1.14 fastqc/0.11.2 nextflow/0.24.1-SNAPSHOT\n";

print CAS "bcl2fastq --barcode-mismatches 0 -o /project/PHG/PHG_Clinical/illumina/$prjid --ignore-missing-positions --no-lane-splitting --ignore-missing-filter --ignore-missing-bcls --runfolder-dir $seqdatadir --sample-sheet /project/PHG/PHG_Clinical/illumina/sample_sheets/$prjid\.bcl2fastq.csv &> /project/PHG/PHG_Clinical/illumina/logs/run_casava_$prjid\.log\n";
print CAS "mkdir /project/PHG/PHG_BarTender/bioinformatics/demultiplexing/$prjid\n" unless (-e "/project/PHG/PHG_BarTender/bioinformatics/demultiplexing/$prjid");
print CAS "cp -R /project/PHG/PHG_Clinical/illumina/$prjid\/Reports /project/PHG/PHG_BarTender/bioinformatics/demultiplexing/$prjid\n" unless (-e "/project/PHG/PHG_BarTender/bioinformatics/demultiplexing/$prjid/Reports");
print CAS "mv /project/PHG/PHG_Clinical/illumina/$prjid\/Stats /project/PHG/PHG_BarTender/bioinformatics/demultiplexing/$prjid\n" unless (-e "/project/PHG/PHG_BarTender/bioinformatics/demultiplexing/$prjid/Stats");

my %completeout; 
my %control;
my %completeout_somatic;

my $prodir = "/project/PHG/PHG_Clinical/processing";
my $outdir = "$prodir\/$prjid/fastq";
my $outnf = "$prodir\/$prjid/analysis";
my $workdir = "$prodir\/$prjid/work";
system("mkdir $prodir\/$prjid") unless (-e "$prodir\/$prjid");
system("mkdir $outdir") unless (-e $outdir);
system("mkdir $outnf") unless (-e $outnf);
system("mkdir $workdir") unless (-e $workdir);
print CAS "cd /project/PHG/PHG_Clinical/processing/$prjid\n";

foreach $dtype (keys %samples) {
  open SSOUT, ">$outdir\/$dtype\.design.txt" or die $!;
  if ($umi) {
    print SSOUT join("\t","SampleID",'SampleID2','SampleName','FamilyID','FqR1','FqR2','BAM','FinalBAM'),"\n";
  }else {
    print SSOUT join("\t","SampleMergeName",'SampleID','SampleName','SubjectID',
		     'FullPathToFqR1','FullPathToFqR2','BAM','OntargetBAM'),"\n";
  }
  my %thash;
  foreach $project (keys %{$samples{$dtype}}) {
    my $datadir =  "/project/PHG/PHG_Clinical/illumina/$prjid/$project/";
    foreach $samp (@{$samples{$dtype}{$project}}) {
      my %info = %{$sampleinfo{$samp}};
      if($info{SubjectID} eq 'GM12878'){ #Positive Control
	  $control{$info{MergeName}}='GM12878';
      }
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
      $thash{$finalrestingplace} = 1;
      $completeout{$info{MergeName}} = $finalrestingplace;
      print CAS "ln -s $datadir/$samp*_R1_*.fastq.gz $finalrestingplace\/$samp\.R1.fastq.gz\n";
      print CAS "ln -s $datadir/$samp*_R2_*.fastq.gz $finalrestingplace\/$samp\.R2.fastq.gz\n";
      print SSOUT join("\t",$info{MergeName},$info{Sample_ID},$info{Sample_Name},
		       $info{SubjectID},"$samp\.R1.fastq.gz","$samp\.R2.fastq.gz",
		       $info{MergeName}.".bam",$info{MergeName}.".final.bam"),"\n";
    }
  }
  foreach my $directory(keys %thash){
    system("mkdir $directory");}
  close SSOUT;
  open TNPAIR, ">$outdir\/$dtype\.design_tumor_normal.txt" or die $!;
  my $tnpairs = 0;
  my $tonlys = 0;
  if ($umi) {
    print TNPAIR join("\t",'TumorID','NormalID','TumorBAM','NormalBAM',
		      'TumorFinalBAM','NormalFinalBAM'),"\n";
  }else {
    print TNPAIR join("\t",'TumorID','NormalID','TumorBAM','NormalBAM',
		      'TumorOntargetBAM','NormalOntargetBAM'),"\n";
  }
  my %sthash;
  foreach my $subjid (keys %samps) {
    my @ctypes = keys %{$samps{$subjid}};
    if ($samps{$subjid}{tumor} && $samps{$subjid}{normal}) {
      print TNPAIR join("\t",$samps{$subjid}{tumor},$samps{$subjid}{normal},
			$samps{$subjid}{tumor}.".bam",
			$samps{$subjid}{normal}.".bam",
			$samps{$subjid}{tumor}.".final.bam",
			$samps{$subjid}{normal}.".final.bam"),"\n";
      $tnpairs ++;
      my $somatic_name=$samps{$subjid}{tumor}."_".$samps{$subjid}{normal};
      my %som_info;
      foreach my $saminfo(keys %sampleinfo){
	if ($saminfo =~ m/$samps{$subjid}{tumor}/){
	  %som_info =%{$sampleinfo{$saminfo}};
	}
      }
      my $finaloutput_somatic = '/project/PHG/PHG_Clinical/'.$som_info{ClinRes};
      unless (-e "$finaloutput_somatic\/$subjid") {
	system("mkdir $finaloutput_somatic\/$subjid");
      }
      my $finalrestingplace_somatic = "$finaloutput_somatic\/$subjid\/$somatic_name";
      if (-e $finalrestingplace_somatic) {
	$finalrestingplace_somatic .= "_".(split(/_|-/,$prjid))[-1];
      }
      $sthash{$finalrestingplace_somatic} = 1;
      $completeout_somatic{$somatic_name} = $finalrestingplace_somatic;
    }
  }
  foreach my $directory (keys %sthash) {
    system("mkdir $directory");
  }
  close TNPAIR;
  my $capture = '/project/shared/bicf_workflow_ref/GRCh38/UTSWV2.bed';
  $capture = '/project/shared/bicf_workflow_ref/GRCh38/MedExome_Plus.bed' if ($dtype eq 'medexomeplus');
  my $mdup = 'picard';
  $mdup = 'fgbio_umi' if ($umi);
  if ($dtype =~ /panel1385|exome|dnaseq/) {
    print CAS "nextflow -C /project/PHG/PHG_Clinical/clinseq_workflows/nextflow.config.super run -w $workdir /project/PHG/PHG_Clinical/clinseq_workflows/alignment.nf --design $outdir\/$dtype\.design.txt --capture $capture --input $outdir --output $outnf --markdups $mdup > $outnf\/$dtype\.nextflow_alignment.log\n";
    print CAS "nextflow -C /project/PHG/PHG_Clinical/clinseq_workflows/nextflow.config.super run -w $workdir /project/PHG/PHG_Clinical/clinseq_workflows/somatic.nf --design $outdir\/$dtype\.design_tumor_normal.txt  --callsvs skip --capture $capture --input $outnf --output $outnf > $outnf\/$dtype\.nextflow_somatic.log &\n" if ($tnpairs);
    print CAS "nextflow -C /project/PHG/PHG_Clinical/clinseq_workflows/nextflow.config.super run -w $workdir /project/PHG/PHG_Clinical/clinseq_workflows/tumoronly.nf --design $outdir\/$dtype\.design.txt --capture $capture --input $outnf --output $outnf > $outnf\/$dtype\.nextflow_tumoronly.log &\n";
    
  }elsif ($dtype =~ m/rnaseq/) {
    print CAS "nextflow -C /project/PHG/PHG_Clinical/clinseq_workflows/nextflow.config.super run -w $workdir /project/PHG/PHG_Clinical/clinseq_workflows/rnaseq.nf --design $outdir\/$dtype\.design.txt --input $outdir --output $outnf --markdups $mdup > $outnf\/$dtype\.nextflow_rnaseq.log\n";
    print CAS "nextflow -C /project/PHG/PHG_Clinical/clinseq_workflows/nextflow.config.super run -w $workdir /project/PHG/PHG_Clinical/clinseq_workflows/tumoronly.nf --design $outdir\/$dtype\.design.txt --genome /project/shared/bicf_workflow_ref/GRCh38/hisat_index --nuctype rna --callsvs skip --capture $capture --input $outnf --output $outnf > $outnf\/$dtype\.nextflow_tumoronly.log &\n";
  }
  print CAS "wait\n";
  unless ($umi) {
    print CAS "cd $outnf\n";
    foreach my $somid(keys %completeout_somatic){
      $finalrestingplace_somatic = $completeout_somatic{$somid};
      print CAS "mv $outnf\/$somid\* $finalrestingplace_somatic\n";
    }
    foreach $sampid (keys %completeout) {
      $finalrestingplace = $completeout{$sampid};
      print CAS "mv $outnf\/$sampid\* $finalrestingplace\n";
    }
    foreach my $posCtrls(keys %control){
      my $prefixName = $posCtrls;
      print CAS "bash /project/PHG/PHG_Clinical/clinseq_workflows/scripts/snsp.sh $prefixName >$prefixName\.snsp\.txt\n";
    }
  }
}
close CAS;
