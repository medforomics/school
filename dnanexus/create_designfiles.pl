#!/usr/bin/perl -w
#run_casava.pl

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my %opt = ();
my $results = GetOptions (\%opt,'help|h','input|i=s','fastqlist|o=s');

open SS, "<$opt{input}" or die $!;
my %rna;
my %dna;
my %panelbed=("heme"=>"UTSW_V4_heme", "pancancer"=>"UTSW_V4_pancancer", "idtrnaseq"=>"UTSW_V4_rnaseq");

while (my $line = <SS>){
  chomp($line);
  $line =~ s/\r//g;
  $line =~ s/ //g;
  $line =~ s/,+$//g;
  if ($line =~ m/^\[Data\]/) {
    $header = <SS>;
    $header =~ s/\r//g;
    chomp($header);
    $header =~ s/Sample_*/Sample_/g;
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
      if ($hash{Sample_Name} =~ m/_Lib/) {
	  $hash{Sample_Name} =~ s/_Lib.*//;
      }
      $hash{Sample_Name} =~ s/T_RNA_panelrnaseq-\d+-\d+/T_RNA_panelrnaseq/;
      $hash{Sample_Project} = $hash{Project} if $hash{Project};
      $hash{Sample_Project} =~ s/\s*$//g;
      $hash{Assay} = lc($hash{Assay});
      unless ($hash{Class}) {
	$hash{Class} = 'tumor';
	$hash{Class} = 'normal' if ($hash{Sample_Name} =~ m/_N_/);
      }
      $hash{SubjectID} = $hash{Sample_Project};
      $hash{Sample_ID} = $hash{Sample_Name};
      my ($ord,$mrn) = split(/-/,$hash{SubjectID});
      if ($hash{Assay} =~ m/rna/) {
	$rna{$ord}{$hash{Sample_Name}} = $panelbed{$hash{Assay}};
      } else {
	$dna{$ord}{lc($hash{Class})}{$hash{Sample_Name}} = $panelbed{$hash{Assay}};
      }
    }
  }
}
my %fastq;
open FQLIST, "<$opt{fastqlist}" or die $!;
while (my $line = <FQLIST>) {
  chomp($line);
  my ($runid,$caseid,$fq)  = split(/\t/,$line);
  my ($sampleid,$rest) = split(/_S\d+_/,$fq);
  #$run{$sampleid} = $runid;
  if ($rest =~ m/R1/) {
    $fastq{$sampleid}{R1} = $fq;
  }else {
    $fastq{$sampleid}{R2} = $fq;
  }
}
open RNA, ">rna_tumor_panel.design.txt" or die $!;
print RNA join("\t","CaseID","SampleID","FqR1","FqR2"),"\n";
foreach my $caseid (keys %rna) {
  foreach my $sid (keys %{$rna{$caseid}}) {
    print RNA join("\t", $caseid,$sid,$fastq{$sid}{R1},$fastq{$sid}{R2}),"\n";
    
  }
}
open DNAS, ">dna_somatic.design.txt" or die $!;
print DNAS join("\t","PanelFile","CaseID","TumorID","TumorFqR1","TumorFqR2","NormalID","NormalFqR1","NormalFqR2"),"\n";
open DNAT, ">dna_tumoronly.design.txt" or die $!;
print DNAT join("\t","PanelFile","CaseID","SampleID","FqR1","FqR2"),"\n";
foreach my $caseid (keys %dna) {
  if ($dna{$caseid}{normal}) {
    my @norms = keys %{$dna{$caseid}{normal}};
    my $nid = shift @norms;
    foreach $sid (keys %{$dna{$caseid}{tumor}}) {
      print DNAS join("\t",$dna{$caseid}{tumor}{$sid},$caseid,$sid,
		      $fastq{$sid}{R1},$fastq{$sid}{R2},
		      $nid,$fastq{$nid}{R1},$fastq{$nid}{R2}),"\n";
    }
    foreach $sid (@norms) {
      print DNAT join("\t",$dna{$caseid}{tumor}{$sid},$caseid,$sid,
		      $fastq{$sid}{R1},$fastq{$sid}{R2}),"\n";
    }
  }else {
    foreach $sid (keys %{$dna{$caseid}{tumor}}) {
      print DNAT join("\t",$dna{$caseid}{tumor}{$sid},$caseid,$sid,
		      $fastq{$sid}{R1},$fastq{$sid}{R2}),"\n";
    }
  }
}
