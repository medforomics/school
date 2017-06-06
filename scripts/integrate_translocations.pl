#!/usr/bin/perl -w
#svvcf2bed.pl

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my %opt = ();
my $results = GetOptions (\%opt,'svcf|v=s','fusion|f=s','help|h');

open OM, "</project/shared/bicf_workflow_ref/GRCh38/utswv2_fusion.keep.txt" or die $!;
while (my $line = <OM>) {
  chomp($line);
  $fusion{$line} = 1;
}
close OM;
open OM, "</project/shared/bicf_workflow_ref/GRCh38/utswv2_known_genefusions.txt" or die $!;
while (my $line = <OM>) {
  chomp($line);
  $known{$line} = 1;
}
close OM;
open OM, "</project/shared/bicf_workflow_ref/GRCh38/genenames.txt" or die $!;
while (my $line = <OM>) {
  chomp($line);
  my ($chr,$start,$end,$ensid,$gene,$type) = split(/\t/,$line);
  $prots{$gene} = 1 if ($type eq 'protein_coding');
}
open OM, "</project/shared/bicf_workflow_ref/GRCh38/utswv2_rearrangement.keep.txt" or die $!;
while (my $line = <OM>) {
  chomp($line);
  $rearrag{$line} = 1;
}
close OM;

if ($opt{fusion} && -e $opt{fusion}) {
  open FUSION, "<$opt{fusion}" or die $!;
  my $header = <FUSION>;
  chomp($header);
  $header =~ s/^#//;
  my @hline = split(/\t/,$header);
  
  while (my $line = <FUSION>) {
    chomp($line);
    my @row = split(/\t/,$line);
    my %hash;
    foreach my $i (0..$#row) {
      $hash{$hline[$i]} = $row[$i];
    }
    my ($left_chr,$left_pos,$left_strand) = split(/:/,$hash{LeftBreakpoint});
    my ($right_chr,$right_pos,$right_strand) = split(/:/,$hash{RightBreakpoint});
    $hash{LeftGene} = (split(/\^/,$hash{LeftGene}))[0];
    $hash{RightGene} = (split(/\^/,$hash{RightGene}))[0];
    next unless ($fusion{$hash{LeftGene}} || $fusion{$hash{RightGene}});
    next unless ($prots{$hash{LeftGene}} && $prots{$hash{RightGene}});
    my $fusionname = join("--",sort {$a cmp $b} ($hash{LeftGene},$hash{RightGene}));
    if ($tloc{$fusionname}) {
      push @{$tloc{$fusionname}{RNAInfo}},[$hash{JunctionReadCount}+$hash{SpanningFragCount},join(":",$left_chr,$left_pos),join(":",$right_chr,$right_pos)];
      $tloc{$fusionname}{SumRNAReads} += $hash{JunctionReadCount}+$hash{SpanningFragCount};
    }else {
      $tloc{$fusionname} = {LeftGene=>$hash{LeftGene},
			    RightGene=>$hash{RightGene},
			    LeftStrand=>$left_strand,
			    RightStrand=>$right_strand,
			    SumRNAReads=>$hash{JunctionReadCount}+$hash{SpanningFragCount},
			    DNAReads=>0};
      push @{$tloc{$fusionname}{RNAInfo}},[$hash{JunctionReadCount}+$hash{SpanningFragCount},join(":",$left_chr,$left_pos),join(":",$right_chr,$right_pos)];
    }
  }
}

my %svs;
my %annot;
open LUMPY, "<$opt{svcf}" or die $!;
while (my $line = <LUMPY>) {
  chomp($line);
    my ($id,$chr,$start,$end,$svtype,$cnv,$numreads,$gene,$exons,$sample,$dgv) = split(/\t/,$line);
  push @{$svs{$id}}, [$chr,$start,$end,$svtype,$numreads,$gene,$exons,$sample];
  $annot{$chr}{$start}{$gene} = $exons;
}

my %svs_gene;
foreach my $id (keys %svs) {
  my @loci = @{$svs{$id}};
  my %gene;
  my %chr;
  my $prevpos = 0;
  my $prevchr = 'chrZ';
  my $toptier = 0;
  my $transloc = 0;
  my $relreads = 0;
  my $allprots = 1;
  my %bnder;
  foreach $locus (sort {$a->[0] cmp $b->[0] && $a->[1] cmp $b->[1]} @loci) {
    my @t1 = @{$locus};
    $relreads ++ if ($rearrag{$t1[5]});
    $allprots = 0 unless $prots{$t1[5]}; 
    if ($t1[3] ne 'BND') {
      next unless ($fusion{$t1[5]});
      push @{$svs_gene{$t1[5]}{$t1[3]}}, [$t1[4],$t1[0],
					 join("-",$t1[1],$t1[2]),$t1[6]];
    }else {
      $bnder{$t1[0]}{$t1[1]}{$t1[5]} = $t1[4];
      $transloc = 1 if ($t1[0] eq $prevchr && $t1[1] - $prevpos > 1e-7);
    }
  }
  next unless ($allprots);
  my @chr = keys %bnder;
  if ($#chr > 0 || $transloc == 1) {
    my @gene;
    my @kloci;
    foreach $chr (keys %bnder) {
      foreach $pos (keys %{$bnder{$chr}}) {
	my @ogen = keys %{$bnder{$chr}{$pos}};
	if (scalar(@ogen) > 1) {
	  if($ogen[0] =~ m/^RP\d+_|^RP\d+-/) {
	    push @kloci, [$chr,$pos,$bnder{$chr}{$pos}{$ogen[1]}];
	    push @gene, $ogen[1];
	  }else {
	    push @kloci, [$chr,$pos,$bnder{$chr}{$pos}{$ogen[0]}];
	    push @gene, $ogen[0]};
	}else {
	  push @kloci, [$chr,$pos,$bnder{$chr}{$pos}{$ogen[0]}];
	  push @gene, $ogen[0];
	}
      }
    }
    my $fusionname = join("--",sort {$a cmp $b} @gene);
    if ($tloc{$fusionname}) {
      $tloc{$fusionname}{DNAReads} += $kloci[0][2];
    }else {
      $tloc{$fusionname} = {LeftGene=>$gene[0],
			    RightGene=>$gene[1],
			    LeftStrand=>'',
			    RightStrand=>'',
			    SumRNAReads=>0,
			    DNAReads=>$kloci[0][2]};
     push @{$tloc{$fusionname}{RNAInfo}},[0,join(":",$kloci[0][0],$kloci[0][1]),join(":",$kloci[1][0],$kloci[1][1])];
    }
  } else {
    next unless $relreads;
    foreach $chr (keys %bnder) {
      foreach $pos (keys %{$bnder{$chr}}) {
	my @genes = keys %{$bnder{$chr}{$pos}};
	next unless $#genes < 1;
	my $gene = $genes[0];
	next unless ($rearrag{$gene});
	push @{$svs_gene{$gene}{'REARG'}}, [$bnder{$chr}{$pos}{$gene},$chr,$pos,$annot{$chr}{$pos}{$gene}];
      }
    }
  }
}

$tloc_out = $opt{svcf};
$tloc_out =~ s/sv.annot.txt/translocations.txt/;
open OUT, ">$tloc_out" or die $!;
print OUT join("\t","FusionName","LeftGene","RightGene","LefttBreakpoint",
	       "RightBreakpoint","LeftStrand","RightStrand","RNAReads",
	       "DNAReads"),"\n";

foreach $fname (keys %tloc) {
  $keepbp = 0;
  if ($known{$fname} && ($tloc{$fname}{SumRNAReads} >= 1 || $tloc{$fname}{DNAReads} >= 3 )) {
    $keepbp = 1;
  }elsif (($tloc{$fname}{SumRNAReads} >= 2 || $tloc{$fname}{DNAReads} >= 20)) {
    $keepbp = 1; 
  }
  next unless ($keepbp);
  foreach $bps (@{$tloc{$fname}{RNAInfo}}) {
    my ($rnareads, $leftbp, $rightbp) = @{$bps};
    print OUT join("\t",$fname,$tloc{$fname}{LeftGene},$tloc{$fname}{RightGene},
		   $leftbp, $rightbp,$tloc{$fname}{LeftStrand},
		   $tloc{$fname}{RightStrand},$rnareads,$tloc{$fname}{DNAReads}),"\n";
  }
}
close OUT;

$reag_out = $opt{svcf};
$reag_out =~ s/sv.annot.txt/svcalls.txt/;
open OUT, ">$reag_out" or die $!;
print OUT join("\t","SVType","Gene","NumberReads","Chr","Loci","Exons"),"\n";
foreach $gene (keys %svs_gene) {
  foreach $svtype (keys %{$svs_gene{$gene}}) {
    my @arrs = sort {$b->[0] <=> $a->[0]} @{$svs_gene{$gene}{$svtype}};
    my ($dnareads,$chr,$pos,$annot) = @{$arrs[0]};
    next unless ( $dnareads > 10);
    print OUT join("\t",$svtype,$gene,$dnareads,$chr,$pos,$annot),"\n";
  }
}
