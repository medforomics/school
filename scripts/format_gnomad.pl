#!/usr/bin/perl -w
#temp.pl

open IN, "</archive/PHG/PHG_Clinical/phg_workflow/analysis/gnomad/gnomad_lcr.final.txt" or die $!;
open OUT, ">gnomad.txt" or die $!;

print OUT join("\t",'CHROM','POS','REF','ALT','GNOMAD_HOM','GNOMAD_AF','AF_POPMAX','GNOMAD_HG19_VARIANT','GNOMAD_LCR'),"\n";
my $header = <IN>;
chomp($header);
my @colnames = split(/\t/,$header);

while (my $line = <IN>) {
  chomp($line);
  my @row = split(/\t/,$line);
  my %hash = (GNOMAD_GEN_HOM=>0,GNOMAD_EX_HOM=>0,GNOMAD_GEN_AN=>0,GNOMAD_EX_AN=>0,GNOMAD_GEN_AC=>0,GNOMAD_EX_AC=>0);
  foreach my $i (0..$#row) {
    $hash{$colnames[$i]} = $row[$i] unless ($row[$i] eq '' || $row[$i] eq '.');
  }
  my $lcr = '0';
  if ($hash{GNOMAD_GEN_LCR} || $hash{GNOMAD_EX_LCR}) {
      $lcr = 1;
  }
  my $hom = 0;
  $hom = $hash{GNOMAD_GEN_HOM} + $hash{GNOMAD_EX_HOM} if ($hash{GNOMAD_GEN_HOM} || $hash{GNOMAD_EX_HOM}); 
  my $an = $hash{GNOMAD_GEN_AN} + $hash{GNOMAD_EX_AN};
  my $ac = $hash{GNOMAD_GEN_AC} + $hash{GNOMAD_EX_AC};
  my @popmaxaf;
  push @popmaxaf, $hash{GNOMAD_GEN_AF_POPMAX} if ($hash{GNOMAD_GEN_AF_POPMAX}); 
  push @popmaxaf, $hash{GNOMAD_EX_AF_POPMAX} if ($hash{GNOMAD_EX_AF_POPMAX}); 
  @popmaxaf = sort {$b <=> $a} @popmaxaf;
  my $popmax = $popmaxaf[0];
  my $af;
  if ($an) {
    $af = sprintf("%.2e",$ac/$an);
  } else {
    next unless $popmax;
    $af = $popmax if ($popmax);
  }
  unless ($popmax) {
    $popmax = $af;
  }
  my $hg19;
  if ($hash{'37_POSITION_GEN'}) {
    $hg19 = $hash{'37_POSITION_GEN'};
  }elsif ($hash{'37_POSITION_EX'}) {
    $hg19 = $hash{'37_POSITION_EX'};
  }
  unless (exists $hash{POSITION} && defined $hg19 && defined $hom && defined $af && defined $popmax) {
    warn "debug\n";
  }
  my ($chr,$pos,$ref,$alt) = split(":",$hash{POSITION});
  unless ($chr =~ m/chr/) {
      $chr = "chr".$chr;
  }
  print OUT join("\t",$chr,$pos,$ref,$alt,$hom,$af,$popmax,$hg19,$lcr),"\n";
}
