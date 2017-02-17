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

open SS, "</project/PHG/PHG_Illumina/sample_sheets/$prjid\.csv" or die $!;
my %sampleinfo;
while (my $line = <SS>){
  chomp($line);
  $line =~ s/\r//g;
  $line =~ s/ //g;
  $line =~ s/,+$//g;
  if ($line =~ m/^\[Data\]/) {
    $header = <SS>;
    $header =~ s/\r//g;
    chomp($header);
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
      $hash{Sample_Project} =~ s/\s*$//g;
      $project = $hash{Sample_Project};
      $samp = $hash{Sample_Name};	
      $dtype = lc($hash{Description});
      $hash{SubjectID} = 'unk' unless $hash{SubjectID};
      $hash{ClinRes} = 'Clinical' unless $hash{ClinRes};
      $sampleinfo{$samp} = \%hash;
      if ($hash{ClinRes} eq 'Clinical') {
	  push @{$samples{$dtype}{$project}}, $samp;
      }
    }
  }
}

open CAS, ">/project/PHG/PHG_Illumina/logs/run_casava_$prjid\.sh" or die $!;
print CAS "#!/bin/bash\n#SBATCH --job-name bcl2fastq\n#SBATCH -N 1\n";
print CAS "#SBATCH -t 1-0:0:00\n#SBATCH -o job_%j.out\n#SBATCH -e job_%j.err\n";
print CAS "#SBATCH --mail-type ALL\n#SBATCH --mail-user brandi.cantarel\@utsouthwestern.edu\n";
print CAS "module load bcl2fastq/2.17.1.14 fastqc/0.11.2 nextflow/0.20.1\n";
print CAS "bcl2fastq --barcode-mismatches 1 --with-failed-reads -o /project/PHG/PHG_Clinical/illumina/$prjid --ignore-missing-positions --min-log-level DEBUG --no-lane-splitting --ignore-missing-bcls --runfolder-dir /project/PHG/PHG_Illumina/BioCenter/$prjid --sample-sheet /project/PHG/PHG_Illumina/sample_sheets/$prjid\.csv &> /project/PHG/PHG_Illumina/logs/run_casava_$prjid\.log\n";
print CAS "/project/PHG/PHG_Illumina/scripts/parse_conversion.stats.pl /project/PHG/PHG_Clinical/illumina/$prjid\n";
foreach $dtype (keys %samples) {
  my $outdir = "/project/PHG/PHG_Clinical\/$dtype\/fastq/";
  my $outnf = "/project/PHG/PHG_Clinical\/$dtype\/analysis";
  open SSOUT, ">$outdir\/samplesheet_$dtype\_$prjid\.txt" or die $!;
  print SSOUT join("\t",'SampleID','SampleName','SubjectID','FullPathToFqR1','FullPathToFqR2'),"\n";
  foreach $project (keys %{$samples{$dtype}}) {
    my $datadir =  "/project/PHG/PHG_Clinical/illumina/$prjid/$project/";
    foreach $samp (@{$samples{$dtype}{$project}}) {
      my %info = %{$sampleinfo{$samp}};
      print CAS "ln -s $datadir/$samp*_R1_*.fastq.gz $outdir\/$samp\.R1.fastq.gz\n";
      print CAS "ln -s $datadir/$samp*_R2_*.fastq.gz $outdir\/$samp\.R2.fastq.gz\n";
      print SSOUT join("\t",$info{Sample_ID},$info{Sample_Name},$info{SubjectID},
		       "$samp\.R1.fastq.gz","$samp\.R2.fastq.gz"),"\n";
    }
  }
  close SSOUT;
  if ($dtype eq 'panel1385') {
    print CAS "nextflow run /project/PHG/PHG_Clinical/clinseq_workflows/germline.nf --design $outdir\/samplesheet_$dtype\_$prjid\.txt --capture /project/shared/bicf_workflow_ref/GRCh38/UTSWV2.bed --input $outdir --output $outnf &> $outnf\/nextflow.log &\n"
  }elsif ($dtype eq 'medexomeplus') {
    print CAS "nextflow run /project/PHG/PHG_Clinical/clinseq_workflows/germline.nf --design $outdir\/samplesheet_$dtype\_$prjid\.txt --capture /project/shared/bicf_workflow_ref/GRCh38/MedExome_Plus.bed --input $outdir --output $outnf &> $outnf\/nextflow.log &\n"
  }elsif ($dtype eq 'panelrnaseq' || $dtype eq 'wholernaseq') {
    print CAS "nextflow run /project/PHG/PHG_Clinical/clinseq_workflows/rnaseq.nf --design $outdir\/samplesheet_$dtype\_$prjid\.txt --input $outdir --output $outnf &> $outnf\/nextflow.log &\n"
  }
}
close CAS;
