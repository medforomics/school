#!/usr/bin/perl 
#migrate_db.pl

my $vcf = shift @ARGV;
my $outfile = $vcf;
$outfile =~ s/int.vcf/uniform.vcf/;
open VCF, "<$vcf" or die $!;
open OUT, ">$outfile" or die $!;

while (my $line = <VCF>) {
  chomp($line);
  if ($line =~ m/#/) {
    print OUT $line,"\n";
    next;
  }
  my ($chrom, $pos,$id,$ref,$alt,$score,
      $filter,$annot,$format,$allele_info) = split(/\t/, $line);
  my %hash = ();
  foreach $a (split(/;/,$annot)) {
    my ($key,$val) = split(/=/,$a);
    $hash{$key} = $val;
  }
  my @deschead = split(/:/,$format);
  my @gtinfo = split(/:/,$allele_info);
  my %gtdata;
  foreach my $i (0..$#deschead) {
    $gtdata{$deschead[$i]} = $gtinfo[$i];
  }
  $newformat = 'GT:DP:AD:AO:RO';
  $gtdata{DP} = $hash{DP} unless ($gtdata{DP});
  unless ($gtdata{AO}) {
    if ($gtdata{AD}){
      ($gtdata{RO},$gtdata{AO}) = split(/,/,$gtdata{AD});
    }elsif ($hash{AD}){
      ($gtdata{RO},$gtdata{AO}) = split(/,/,$hash{AD});
    }elsif ($gtdata{NR} && $gtdata{NV}) {
      $gtdata{DP} = $gtdata{NR}; 	
      $gtdata{AO} = $gtdata{NV};
      $gtdata{RO} = $gtdata{DP} - $gtdata{AO};
    }elsif ($hash{TR}) {
      $gtdata{AO} = $hash{TR};
      $gtdata{DP} = $hash{TC};
      $gtdata{RO} = $gtdata{DP} - $gtdata{AO};
    }
  }
  if ($gtdata{DP} < $gtdata{AO}+$gtdata{RO}) {
      warn "inconsistent\n";
  }
  unless ($gtdata{AD}) {
    if (exists $gtdata{RO} && exists $gtdata{AO}) {
      $gtdata{AD} = join(",",$gtdata{RO},$gtdata{AO});
    }
  }
  unless (exists $gtdata{GT} && exists $gtdata{DP} && exists $gtdata{AD} &&  exists $gtdata{AO} && exists $gtdata{RO}) {
    warn "Missing Information\n";
  }
  $newgt = join(":",$gtdata{GT},$gtdata{DP},$gtdata{AD},$gtdata{AO},$gtdata{RO});
  print OUT join("\t",$chrom,$pos,$id,$ref,$alt,$score,$filter,$annot,$newformat,$newgt),"\n";
}
