#!/usr/bin/perl -w
#run_casava.pl

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my %opt = ();
my $results = GetOptions (\%opt,'help|h','input|i=s','output|o=s');

open SS, "<$opt{input}" or die $!;
open SSOUT, ">$opt{output}" or die $!;

my %spairs;
my %stype;
my %samples;

while (my $line = <SS>){
  chomp($line);
  $line =~ s/\r//g;
  $line =~ s/ //g;
  $line =~ s/,+$//g;
  if ($line =~ m/^\[Data\]/) {
    print SSOUT join("\n","[Settings]","ReverseComplement,0","Read2UMILength,8","TrimUMI,1"),"\n";
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
