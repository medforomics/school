#!/usr/bin/perl -w
#cosmicbkpnt2bed.pl

#load bedtool
my $refidx = shift @ARGV;
open IN, "gunzip -c CosmicBreakpointsExport.tsv.gz |" or die $!;
open OUT, ">cosmic_sv.unsort.bed" or die $!;

my $header = <IN>;
chomp($header);
$header =~ s/ //g;
@colnames = split(/\t/,$header);

while ($line = <IN>) {
  chomp($line);
  my @row = split(/\t/, $line);
  my %hash;
  foreach my $i (0..$#row) {
    $hash{$colnames[$i]} = $row[$i];
  }
  next if ($hash{ChromTo} > 23 || $hash{ChromFrom} > 23);
  $hash{ChromFrom} = 'chr'.$hash{ChromFrom};
  $hash{ChromTo} = 'chr'.$hash{ChromTo};
  if ($hash{ChromFrom} eq 'chr23') {
    $hash{ChromFrom} = 'chrX';
  }if ($hash{ChromTo} eq 'chr23') {
    $hash{ChromTo} = 'chrX';
  }
  my $split = 1;
  if ($hash{LocationTomax} < $hash{LocationTomin}) {
    $temp = $hash{LocationTomax};
    $hash{LocationTomax} = $hash{LocationTomin};
    $hash{LocationTomin} = $temp;
  }
  if ($hash{LocationFrommax} < $hash{LocationFrommin}) {
    $temp = $hash{LocationFrommax};
    $hash{LocationFrommax} = $hash{LocationFrommin};
    $hash{LocationFrommin} = $temp;
  }
  if ($hash{LocationFrommax} < $hash{LocationFrommin}) {
    $temp = $hash{LocationFrommax};
    $hash{LocationFrommax} = $hash{LocationFrommin};
    $hash{LocationFrommin} = $temp;
  }
  if ($hash{MutationType} =~ m/deletion|insertion|intrachromosomal amplicon|inverted/) {
    if ($hash{LocationTomax} - $hash{LocationFrommin} <= 10000) {
      $split = 0;
    }
  }
  if ($split) {
    unless ($done{$hash{ChromFrom}}{$hash{LocationFrommin}}{$hash{LocationFrommax}}) {
      if ($hash{LocationFrommin} == $hash{LocationFrommax}) {
	$hash{LocationFrommin} --;
      }
      print OUT join("\t",$hash{ChromFrom},$hash{LocationFrommin},$hash{LocationFrommax}),"\n";
      $done{$hash{ChromFrom}}{$hash{LocationFrommin}}{$hash{LocationFrommax}} = 1;
    }
    unless ($done{$hash{ChromTo}}{$hash{LocationTomin}}{$hash{LocationTomax}}) {
      if ($hash{LocationTomin} == $hash{LocationTomax}) {
	$hash{LocationTomin} --;
      }
      print OUT join("\t",$hash{ChromTo},$hash{LocationTomin},$hash{LocationTomax}),"\n";
      $done{$hash{ChromTo}}{$hash{LocationTomin}}{$hash{LocationTomax}} = 1;
    }
  } else {
    next if ($done{$hash{ChromFrom}}{$hash{LocationFrommin}}{$hash{LocationTomax}});
    if ($hash{LocationTomax} < $hash{LocationFrommin}) {
      $temp = $hash{LocationTomax};
      $hash{LocationTomax} = $hash{LocationFrommin};
      $hash{LocationFrommin} = $temp;
    }
    print OUT join("\t",$hash{ChromFrom},$hash{LocationFrommin},$hash{LocationTomax}),"\n";
    $done{$hash{ChromFrom}}{$hash{LocationFrommin}}{$hash{LocationTomax}} = 1;
  }
}
system(qq{sort -V -k 1,1 -k 2,2n cosmic_sv.unsort.bed > cosmic_sv.unmerged.bed});
system(qq{bedtools merge -i cosmic_sv.unmerged.bed | bedtools sort -faidx $refidx > cosmic_sv.bed});
