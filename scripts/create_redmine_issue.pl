#!/usr/bin/perl -w
#samplesheet2ticket.pl

use JSON;
use Redmine::API;

my $in = shift @ARGV;
open SS, "<$in" or die $!;

my ($author, $subject, $created);
my %samples;

while (my $line = <SS>){
  chomp($line);
  $line =~ s/\r//g;
  $line =~ s/,+$//g;
  if ($line =~ m/Experiment/) {
    ($field,$subject) = split(/,/,$line);
  }
  if ($line =~ m/Investigator/) {
    ($field,$author) = split(/,/,$line);
  }
  if ($line =~ m/^Date/) {
    ($field,$created) = split(/,/,$line);
    my ($year,$month,$day) = split(/-/,$created);
  }
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
      $hash{SampleName} =~ s/_Lib.*//g;
      my ($orderid,$mrn) = split(/-/,$hash{Project});
      push @{$samples{$hash{Project}}}, $hash{SampleName};
    }
  }
}

my $descr = join("\n",keys %samples);

my %project = (id=>125,'name'=>'CLIA Lab Orders');
my %status = (id=>1,'name'=>'New');
my %tracker = ('id'=>18,'name'=>'CLIATask');
my @list = ("Check NuCLIA for PASS QC","Check pipeline.txt for FAIL",
	    "Confirm which samples have normal and RNA Samples");
my @checklist;
foreach my $l (@list) {
  my %lhash = ("is_done"=>0,"subject"=>$l);
  push @checklist, \%lhash;
}
my %ticket_info=(project_id=>125,status_id=>1,tracker_id=>18,
		 subject=>$subject,description=>$descr,assigned_to_id=>55,
		 "checklists_attributes"=>\@checklist);

my $c = Redmine::API->new('auth_key' => '044bc8137f9e05ce9a2cb9f300099aab41ef667e', base_url => 'bicf.redmine.swmed.edu', trace => 1);

$run_issue = $c->issues->issue->create(%ticket_info);
my %runticket = %{$run_issue};

foreach $prjid (keys %samples) {
  my $tdescr = join("\n",@{$samples{$prjid}});
  my %case=(project_id=>125,status_id=>1,tracker_id=>18,
	    parent_issue_id=>$runticket{issue}{id},subject=>$prjid,
	    description=>$tdescr,assigned_to_id=>55,
	   );
  $case_issue = $c->issues->issue->create(%case);
  print $case_issue->{issue}{id},"\n";
}
