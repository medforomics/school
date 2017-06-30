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
      system("perl /project/PHG/PHG_Clinical/clinseq_workflows/scripts/init_workflows.pl -p $prjid");
      system("source /etc/profile.d/modules.sh");
      system("sbatch /project/PHG/PHG_Clinical/illumina/logs/run_casava_$prjid\.sh");
  }
}
