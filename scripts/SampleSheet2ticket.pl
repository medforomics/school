#!/usr/bin/perl -w
#samplesheet2ticket.pl

my $in = shift @ARGV;

my ($prefix,@ext) = split(/\./, $in);
my $out = $prefix.".ticket.csv";

open SS, "<$in" or die $!;
open OUT, ">$out" or die $!;
print OUT join(",","Tracker","Status","Priority","Subject","Author","Description",
	       "Category","Start date","Parent"),"\n";

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
      push @{$samples{$hash{Project}}}, $hash{SampleName};
    }
  }
}

my $descr = join(" ",keys %samples);

print OUT join(",","CLIA Lab Orders","New","Normal",$subject,$author,
	       $descr,'CompPanCancer',$created,''),"\n";

foreach $prjid (keys %samples) {
  my $tdescr = join(" ",@{$samples{$prjid}});
  print OUT join(",","CLIA Lab Orders","New","Normal",$prjid,$author,
		 $tdescr,'CompPanCancer',$created,1),"\n";
}
