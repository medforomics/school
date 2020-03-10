#!/usr/bin/perl


use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %opt = ();
my $results = GetOptions (\%opt,'in|i=s','out|o=s','genelist|g=s','tid|t=s');

unless  ($opt{in} && $opt{genelist}) {
    print "vcf2bed.pl <vcf_file> <outfile>\n";
    exit 1;
}

open OM, "<$opt{genelist}" or die $!;
while (my $line = <OM>) {
  chomp($line);
  $keep{$line} = 1;
}
close OM;

my $outfile;
if ($opt{out}) {
    $outfile = $opt{out};
}else {
    $outfile = (split(/.vcf/,$opt{in}))[0].".bed";
}

my $vcfFile = $opt{in};
open VCF, "gunzip -c $vcfFile | "  or die $!;
open OUT, ">$outfile" or die $!;

my @gtheader;
W1:while (my $line = <VCF>) {
  chomp($line);
  if ($line =~ m/#/) {
    if ($line =~ m/#CHROM/) {
      ($chrom, $pos,$id,$ref,$alt,$score,
       $filter,$info,$format,@gtheader) = split(/\t/, $line);
      print OUT join("\t",'#CHROM','From','To','Filter',@gtheader),"\n";
    }
    next;
  }
  my ($chrom, $pos,$id,$ref,$alt,$score,
      $filter,$annot,$format,@gts) = split(/\t/, $line);
  my %hash = ();
  foreach $a (split(/;/,$annot)) {
    my ($key,$val) = split(/=/,$a);
    $hash{$key} = $val;
  }
  next if ($hash{CHR2} && $hash{CHR2} ne $chrom);
  my $startPos = $pos;
  unless ($hash{'END'}) {
    next unless ($alt =~ m/^[ATGC]+$/ && $ref =~ m/^[ATGC]+$/);
    if (length($ref) > length($alt)) {
      my $diff = substr($ref, length($alt));
      $hash{'END'} = $startPos + length($diff);
      $hash{SVTYPE} = 'DEL';

    } elsif (length($alt) > length($ref)) {
      $startPos --;
      $hash{'END'} = $startPos + length($alt);
      $hash{SVTYPE} = 'INS';
    }
  }
  @deschead = split(/:/,$format);
  foreach my $i (0..$#gtheader) {
      if ($gtheader[$i] eq $opt{tid}) {
	  my @gtinfo = split(/:/,$gts[$i]);
	  my %gtdata;
	  foreach my $i (0..$#deschead) {
	      $gtdata{$deschead[$i]} = $gtinfo[$i];
	  }
	  next W1 unless ($gtdata{AO} && $gtdata{AO} > 20);
      }
  }
  next if ($hash{SVTYPE} eq 'BND');
  my %glist;
  foreach $trx (split(/,/,$hash{ANN})) {
    my ($allele,$effect,$impact,$gene,$geneid,$feature,
	$featureid,$biotype,$rank,$codon,$aa,$cdna_pos,$len_cdna,
	$aapos,$distance,$err) = split(/\|/,$trx);
    next unless $gene=~ m/\w+/;
    $glist{$gene} = 1 if $keep{$gene};
  }
  next unless (scalar(keys %glist) > 1);
  print OUT join("\t",$chrom,$startPos,$hash{'END'},join(":",$chrom,$pos,$id,$hash{SVTYPE}),$filter,join(";",keys %glist),$format,@gts),"\n";
}

close VCF;
