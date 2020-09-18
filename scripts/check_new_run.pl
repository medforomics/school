#!/usr/bin/perl 
#check_new_run.pl

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my %opt = ();
my $results = GetOptions (\%opt,'help|h','split|i=s','test|o=s',
			  'xmlcheck|u=s','cloud|c');
my $opts;
if ($opt{split}) {
    $opts .= " -s 1"
}if ($opt{test}) {
    $opts .= " -t 1"
}if ($opt{xmlcheck}) {
    $opts .= " -x 1"
}
my $refdir = "/project/shared/bicf_workflow_ref/human/grch38_cloud/dnaref";
my @scriptpath = split(/\//,$0);
pop @scriptpath;
pop @scriptpath;
my $execdir = join("/",@scriptpath);

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
      if ($opt{cloud}) {
	  system(qq{cp $execdir\/dnanexus/init_dnanexus.sh $procbase});
	  open OUT, ">$procbase\/run_$prjid.sh" or die $!;
	  print OUT qq{#!/bin/bash\n$procbase\/init_dnanexus.sh -p $prjid -b $execdir -c $prodir $opts\n};
	  close OUT;
	  system("sbatch -p 32GB $procbase\/run_$prjid.sh");
      } else {
	  system(qq{$execdir\/utsw_biohpc/create_directories.sh -p $prjid -c $prodir});
	  system(qq{cp $execdir\/utsw_biohpc/init_workflows.sh $procbase});
	  open OUT, ">$procbase\/run_$prjid.sh" or die $!;
	  print OUT qq{#!/bin/bash\n$procbase\/init_workflows.sh -p $prjid -b $execdir -c $prodir -r $refdir $opts\n};
	  close OUT;
	  system("sbatch -p 32GB $procbase\/run_$prjid.sh");
      }
  }
}
