#!/usr/bin/perl 
#check_new_run.pl

my @samplesheets = `/bin/ls /project/PHG/PHG_BarTender/bioinformatics/sample_sheets/*.csv`;
chomp(@samplesheets);
foreach my $ss (@samplesheets) {
  my @filepath = split(/\//,$ss);
  my $fname = pop @filepath;
  my $prjid = (split(/\./,$fname))[0];
  unless (-e "/project/PHG/PHG_Clinical/illumina/sample_sheets/$fname") {
      system("cp $ss /project/PHG/PHG_Clinical/illumina/sample_sheets/$fname");
      system("source /etc/profile.d/modules.sh");
      open OUT, ">run_$prjid.sh" or die $!;
      print OUT qq{#!/bin/bash\n/project/PHG/PHG_Clinical/devel/clinseq_workflows/scripts/init_workflows.sh -p $prjid -r /project/shared/bicf_workflow_ref/human/GRCh38\n};
      close OUT;
      system("sbatch -p 32GB run_$prjid.sh");
  }
}
