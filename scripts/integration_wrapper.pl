#!/usr/bin/perl -w
#integration_wrapper.pl
#use Time::localtime;

my ($month,$year) = (localtime)[4,5];
$year += 1900;
$month++;
$date=$year.sprintf("%02s",$month);

my ($refdir) = $ARGV[0];
my @splitPath = split(/\//, $refdir);
my $subject =$splitPath[-1];
my @directories = `find $refdir -type d| sed 's|$refdir/||'`;

my ($fusion,$svcf)=("") x 2;
my $tumorid="no_tumor";
my $normalid="no_normal";
my $somaticid="no_somatic";
my $rnaseqid="no_rnaseq";

foreach my $directory(@directories){
  chomp $directory;
  if($directory =~ m/-T_DNA_panel1385/){
    if($directory !~ m/-N_DNA_panel1385/){
      $tumorid=$directory;
      $svcf=$refdir."/".$tumorid."/".$tumorid.".sv.annot.txt";}}
  elsif($directory =~ m/-N_DNA_panel1385/){
    if($cirectory !~ m/-T_DNA_panel1385/){
      $normalid=$directory;}}
  elsif($directory =~ m/_RNA_panelrnaseq/){
      $rnaseqid=$directory;
      $fusion=$refdir."/".$rnaseqid."/".$rnaseqid.".starfusion.txt";}
}
if($tumorid ne "no_tumor" && $normalid ne "no_normal"){
  $somaticid =$tumorid."_".$normalid;}

system("perl /project/PHG/PHG_Clinical/clinseq_workflows/scripts/integrate_vcfs.pl $subject $subject $tumorid $somaticid $rnaseqid");
system("perl /project/PHG/PHG_Clinical/clinseq_workflows/scripts/integrate_translocations.pl -svcf $svcf -fusion $fusion");
system("python /project/PHG/PHG_Clinical/clinseq_workflow/IntellispaceDemographics/gatherdemographics.py -i $subject -u phg_workflow -p UGMP_Cl1nS3q -o /project/PHG/PHG_Sap/input/$date/$subject/$subject.xml");
system("ln -s /project/shared/bicf_workflow_ref/vcf2maf/.vep ~/.vep");
system("zcat $subject.philips.vcf.gz | java -jar /cm/shared/apps/snpeff/4.2/SnpSift.jar filter \"(FILTER = 'PASS')\"  > $subject.pass.vcf");
system("perl /project/shared/bicf_workflow_ref/vcf2maf/vcf2maf.pl --input $subject.pass.vcf --output $subject.maf --species homo_sapiens --ncbi-build GRCh38 --ref-fasta /project/shared/bicf_workflow_ref/vcf2maf/.vep/homo_sapiens/87_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa --cache-version 87 --vep-path /project/shared/bicf_workflow_ref/vcf2maf/variant_effect_predictor --tumor-id $tumorid");




