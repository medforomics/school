#!/usr/bin/perl -w
#run_casava.pl

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my %opt = ();
my $results = GetOptions (\%opt,'help|h','prjid|p=s');

if (!defined $opt{prjid} || $opt{help}) {
  $usage = <<EOF;
usage: $0 -p prjid

-p prjid -- this is the project name in /project/PHG/PHG_Illumina/BioCenter/ 140505_SN7001189_0117_AH7LRLADXX

EOF
  die $usage,"\n";
}

my @execdir = split(/\//,$0);
pop @execdir;
pop @execdir;
$baseDir = join("/",@execdir);

my $prjid = $opt{prjid};
my $oriss = "/project/PHG/PHG_Clinical/illumina/sample_sheets/$prjid\.csv";
my $newss = "/project/PHG/PHG_Clinical/illumina/sample_sheets/$prjid\.bcl2fastq.csv";
my $capturedir = "/project/shared/bicf_workflow_ref/GRCh38/clinseq_prj";

my $seqdatadir = "/project/PHG/PHG_Illumina/BioCenter/$prjid";
if (-e "/project/PHG/PHG_Illumina/Research/$prjid") {
  $seqdatadir = "/project/PHG/PHG_Illumina/Research/$prjid";
}

$umi = `grep "<Read Number=\\\"2\\\" NumCycles=\\\"14\\\" IsIndexedRead=\\\"Y\\\" />" $seqdatadir/RunInfo.xml`;

open SS, "<$oriss" or die $!;
open SSOUT, ">$newss" or die $!;

my %sampleinfo;
my %stype;
my %spairs;
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
      $hash{Sample_Project} = $hash{Project} if $hash{Project};
      $hash{Sample_Project} =~ s/\s*$//g;
      $hash{Assay} = 'panel1385' if ($hash{Assay} eq 'dnaseqdevelopment');
      $hash{Assay} = 'panel1385' if ($hash{Assay} eq 'panel1385');
      $hash{Assay} = 'panel1385v2' if ($hash{MergeName} =~ m/panel1385v2/);
      $hash{Assay} = 'panelrnaseq' if ($hash{MergeName} =~ m/panelrnaseq/);
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
      my $clinres = 'complete';
      if (($hash{Description} && $hash{Description} =~ m/research/i) ||
	  ($hash{Sample_Name} !~ m/ORD/ && $hash{SubjectID} !~ m/GM12878|ROS/)) {
	$clinres = 'toresearch';
      }
      $hash{ClinRes} = $clinres;
      $hash{Sample_ID} = $hash{Sample_Name};
      $stype{$hash{SubjectID}} = $clinres;
      $spairs{$hash{SubjectID}}{lc($hash{Class})}{$hash{MergeName}} = 1;
      $sampleinfo{$hash{Sample_Name}} = \%hash;
      push @{$samples{lc($hash{Assay})}{$hash{SubjectID}}}, $hash{Sample_Name};

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
print CAS "module load bcl2fastq/2.17.1.14 fastqc/0.11.2 nextflow/0.27.6\n";

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

open TNPAIR, ">$outdir\/design_tumor_normal.txt" or die $!;
my $tnpairs = 0;
print TNPAIR join("\t",'PairID','TumorID','NormalID','TumorBAM','NormalBAM',
		  'TumorFinalBAM','NormalFinalBAM'),"\n";
foreach my $subjid (keys %spairs) {
  my @ctypes = keys %{$spairs{$subjid}};
  if ($spairs{$subjid}{tumor} && $spairs{$subjid}{normal}) {
    my @tumors = keys %{$spairs{$subjid}{tumor}};
    my @norms = keys %{$spairs{$subjid}{normal}};
    my $pct = 0;
    foreach $tid (@tumors) {
      foreach $nid (@norms) {
	my $pair_id = $subjid;
	if ($pct > 1) {
	  $pair_id .= ".$pct";
	}
	print TNPAIR join("\t",$pair_id,$tid,$nid,$tid.".bam",$nid.".bam",
			  $tid.".final.bam",$nid.".final.bam"),"\n";
	$pct ++;
	$tnpairs ++;
      }
    }
  }
}
close TNPAIR;

foreach $dtype (keys %samples) {
  open SSOUT, ">$outdir\/$dtype\.design.txt" or die $!;
  print SSOUT join("\t","SampleID",'SampleID2','SampleName','FamilyID','FqR1','FqR2','BAM','FinalBAM'),"\n";
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
      unless (-e $finalrestingplace) {
	system("mkdir $finalrestingplace");
      }
      print CAS "ln -s $datadir/$samp*_R1_*.fastq.gz $finalrestingplace\/$samp\.R1.fastq.gz\n";
      print CAS "ln -s $datadir/$samp*_R2_*.fastq.gz $finalrestingplace\/$samp\.R2.fastq.gz\n";
      print SSOUT join("\t",$info{MergeName},$info{Sample_ID},$info{Sample_Name},
		       $info{SubjectID},"$samp\.R1.fastq.gz","$samp\.R2.fastq.gz",
		       $info{MergeName}.".bam",$info{MergeName}.".final.bam"),"\n";
    }
  }
  close SSOUT;
  my $capture = "$capturedir\/UTSWV2.bed";
  $capture = "$capturedir\/MedExome_Plus.bed" if ($dtype eq 'medexomeplus');
  $capture = "$capturedir\/UTSWV2_2.bed" if ($dtype eq 'panel1385v2');
  my $mdup = 'picard';
  $mdup = 'fgbio_umi' if ($umi);
  $mdup = 'skip' if ($dtype =~ m/panelrnaseq/);
  my $germopts = '';
  if ($dtype =~ /panel1385|exome|dnaseq/) {
    my $alignwf = "$baseDir\/alignment.nf";
    unless ($umi) {
      $alignwf = "$baseDir\/alignmentV1.nf";
    }
    print CAS "nextflow -C $baseDir\/nextflow.config run -w $workdir $alignwf --design $outdir\/$dtype\.design.txt --capture $capture --input $outdir --output $outnf --markdups $mdup > $outnf\/$dtype\.nextflow_alignment.log\n";
  } elsif ($dtype =~ m/rnaseq/) {
    print CAS "nextflow -C $baseDir\/nextflow.config run -w $workdir $baseDir\/rnaseq.nf --design $outdir\/$dtype\.design.txt --input $outdir --output $outnf --markdups $mdup > $outnf\/$dtype\.nextflow_rnaseq.log\n";
    $germopts = " --genome /project/shared/bicf_workflow_ref/GRCh38/hisat_index --nuctype rna --callsvs skip"
  }
  foreach $project (keys %spairs) {
      foreach $class (keys  %{$spairs{$project}}) {
	  foreach $samp (keys %{$spairs{$project}{$class}}) {
	      print CAS "mv $outnf\/$samp\.* $outnf\/$samp\_* $outnf\/$project\/$samp\n";
	  }
      }
  }
  print CAS "ln -s $outnf\/*/*/*.bam $outnf\n";  ####check me out
  print CAS "nextflow -C $baseDir\/nextflow.config run -w $workdir $baseDir\/tumoronly.nf --design $outdir\/$dtype\.design.txt $germopts --capture $capture --input $outnf --output $outnf > $outnf\/$dtype\.nextflow_tumoronly.log\n";
}
print CAS "nextflow -C $baseDir\/nextflow.config run -w $workdir $baseDir\/somatic.nf --design $outdir\/design_tumor_normal.txt  --callsvs skip --input $outnf --output $outnf > $outnf\/nextflow_somatic.log &\n" if ($tnpairs);
print CAS "wait\n";
print CAS "cd $outnf\n";

foreach $case (keys %stype) {
    print CAS "rsync -avz $case /project/PHG/PHG_Clinical/".$stype{$case},"\n";
  print CAS "rsync -avz --exclude=\"*bam*\" $case /project/PHG/PHG_BarTender/bioinformatics/seqanalysis/".$stype{$case},"\n" if ($stype{$case} eq 'complete');
}

foreach $project (keys %spairs) {
  foreach $class (keys  %{$spairs{$project}}) {
    foreach $samp (keys %{$spairs{$project}{$class}}) {
      print CAS "curl \"http://nuclia-test.biohpc.swmed.edu:8080/NuCLIAVault/addPipelineResults?token=\$nucliatoken&subjectName=$project&sampleName=$samp&runName=$opt{prjid}\"\n";
    }
  }
}
foreach my $posCtrls(keys %control){
  my $prefixName = $posCtrls;
  print CAS "bash $baseDir\/scripts/snsp.sh $prefixName >$prefixName\.snsp\.txt\n";
}
close CAS;
