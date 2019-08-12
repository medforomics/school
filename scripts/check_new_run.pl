#!/usr/bin/perl 
#check_new_run.pl

my $refdir = "/project/shared/bicf_workflow_ref/human/GRCh38";
my @scriptpath = split(/\//,$0);
pop @scriptpath;
my $execdir = join("/",@scriptpath);

exit;

my @samplesheets = `/bin/ls /project/PHG/PHG_BarTender/bioinformatics/sample_sheets/*.csv`;
chomp(@samplesheets);
foreach my $ss (@samplesheets) {
  my @filepath = split(/\//,$ss);
  my $fname = pop @filepath;
  my $prjid = (split(/\./,$fname))[0];
  my $prodir = "/project/PHG/PHG_Clinical/processing";
  my $procbase = "$prodir\/$prjid";
  unless (-e "/project/PHG/PHG_Clinical/illumina/sample_sheets/$fname") {
      system("cp $ss /project/PHG/PHG_Clinical/illumina/sample_sheets/$fname");
      system("source /etc/profile.d/modules.sh");
      system(qq{$execdir\/create_directories.sh -p $prjid -c $prodir});
      system(qq{cp $execdir\/init_workflows.sh $procbase});
      open OUT, ">run_$prjid.sh" or die $!;
      print OUT qq{#!/bin/bash\n$procbase\/init_workflows.sh -p $prjid -b $execdir -c $prodir -r $refdir -x\n};
      close OUT;
      system("sbatch -p 32GB run_$prjid.sh");
  }
}
