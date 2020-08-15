#!/usr/bin/perl -w
#validation_json2txt.pl

use JSON;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my %opt = ();
my $results = GetOptions (\%opt,'analdir|i=s','genelist|l=s');

$analdir="cases";
if ($opt{analdir}) {
    $analdir=$opt{analdir};
}
$genelist='/project/shared/bicf_workflow_ref/human/grch38_cloud/panels/UTSW_V4_heme/genelist.txt';
if ($opt{genelist}) {
    $genelist=$opt{genelist};
}

my %keep;
open OM, "<$genelist" or die $!;
while (my $line = <OM>) {
    chomp($line);
    $keep{$line} = 1;
}

my @vtypes = ('variants','cnvs','translocations');
my %cnvs;
my %fusions;
my %vars;
my %cases;
foreach $jsonfile (@ARGV) {
  $jsontxt = `cat $jsonfile`;
  $jsonref  = decode_json($jsontxt);
  $prefix = (split(/\.|\//,$jsonfile))[1];
  $cases{$prefix} = 1;
  foreach my $vtype (@vtypes) {
    my %csam;
    foreach $tref (@{$jsonref->{$vtype}}) {
      my %hash = %{$tref};
      if ($vtype eq 'variants') {
	next unless $keep{$hash{geneName}};
	$vars{$prefix}{$hash{chrom}}{$hash{pos}} = \%hash;
      }elsif ($vtype eq 'cnvs') {
	my @genes = @{$hash{genes}};
	foreach my $gid (@genes) {
	  next unless $keep{$gid};
	  $cnvs{$prefix}{$gid} = \%hash;
	}
      }elsif ($vtype eq 'translocations') {
	$fusions{$prefix}{$hash{leftGene}}{$hash{rightGene}} = %hash;
      }
    }
  }
}

open SNV, ">snvindel_validation.txt" or die $!;
print SNV join("\t",'Case','ID','CHR','POS','GENE','AA','E-TAF','E-NAF',
		'O-TAF','TDP','O-NAF','NDP','DIF-AF'),"\n";
open FUS, ">fusion_validation.txt" or die $!;
print FUS join("\t",'CaseID','Reported-FusionName',
	       'Reported-RNA_CT',@header),"\n";

foreach $caseID (keys %cases) {
  my $fusionfile = $analdir."/".$caseID."/".$caseID.".translocations.answer.txt";
  open IN, "<$fusionfile" or die $!;
  my $head = <IN>;
  chomp($head);
  my @colnames = split(/\t/,$head);
  my %results;
  while (my $line = <IN>) {
    chomp($line);
    my @row = split(/\t/,$line);
    my %hash;
    foreach my $i (0..$#row) {
      $hash{$colnames[$i]} = $row[$i];
    }
    $results{$hash{LeftGene}}{$hash{RightGene}} = 1;
    if ($fusions{$caseID}{$hash{LeftGene}}{$hash{RightGene}}) {
      print FUS join("\t",$caseID,$fusions{$caseID}{$hash{LeftGene}}{$hash{RightGene}}{fusionName},$fusions{$caseID}{$hash{LeftGene}}{$hash{RightGene}}{rnaReads},$line),"\n";
    }
  }
  foreach $lg (keys %{$fusions{$caseID}}) {
    foreach $rg (keys %{$fusions{$caseID}{$lg}}) {
      next if $results{$lg}{$rg};
      print FUS join("\t",$caseID,$fusions{$caseID}{$lg}{$rg}{fusionName},
		     $fusions{$caseID}{$lg}{$rg}{rnaReads}),"\n";
    }
  }
  close IN;
  my $file = `ls $analdir/$caseID/$caseID.set1.vcf.gz`;
  open IN, "gunzip -c $file |" or die $!;
  my $tumorid;
  my $normalid;
  my %found;
 W1:while (my $line = <IN>) {
    chomp($line);
    if ($line =~ m/^#CHROM/) {
      my @header = split(/\t/,$line);
      ($chrom, $pos,$id,$ref,$alt,$score,
       $filter,$info,$format,@gtheader) = split(/\t/, $line);
      $tumorid = $gtheader[0];
      ($normalid,@o) = grep(/N_DNA/,@gtheader);
      @sampids = @gtheader;
      next;
    }
    next if ($line =~ m/^#/);
    my ($chrom, $pos,$id,$ref,$alt,$score,
	$filter,$annot,$format,@gts) = split(/\t/, $line);
    next if ($ref =~ m/\./ || $alt =~ m/\./ || $alt=~ m/,X/);
    next unless $vars{$caseID}{$chrom}{$pos};
    my %hash = ();
    foreach $a (split(/;/,$annot)) {
      my ($key,$val) = split(/=/,$a);
      $hash{$key} = $val unless ($hash{$key});
    }
    my %gtinfo = ();
    my @deschead = split(/:/,$format);
  F1:foreach my $k (0..$#gtheader) {
      my $subjid = $gtheader[$k];
      my $allele_info = $gts[$k];
      my @ainfo = split(/:/, $allele_info);
      my @mutallfreq = ();
      foreach my $k (0..$#ainfo) {
	$gtinfo{$subjid}{$deschead[$k]} = $ainfo[$k];
      }
      next F1 unless ($gtinfo{$subjid}{DP} && $gtinfo{$subjid}{DP} ne '.'
		      && $gtinfo{$subjid}{DP} >= 1);
      my @altct = split(/,/,$gtinfo{$subjid}{AD});
      my $refct = shift @altct;
      @altct2 = split(/,/,$gtinfo{$subjid}{AO});
      my $total = $refct;
      foreach  my $act (@altct) {
	next if ($act eq '.');
	$total += $act;
	push @mutallfreq, sprintf("%.4f",$act/$gtinfo{$subjid}{DP});
      }
      $gtinfo{$subjid}{MAF} = $mutallfreq[0];
    }
    my %val = %{$vars{$caseID}{$chrom}{$pos}};
    $found{$chrom}{$pos} = 1;
    my $tmaf = $gtinfo{$tumorid}{MAF};
    $val{'normalAltFrequency'} = 0 unless (exists $val{'normalAltFrequency'});
    print SNV join("\t",$caseID,$chrom,$pos,$val{'geneName'},$val{'notation'},
		   sprintf("%.2f",$val{'tumorAltFrequency'}),
		   sprintf("%.2f",$val{'normalAltFrequency'}),
		   $tmaf,$gtinfo{$tumorid}{DP},
		   $gtinfo{$normalid}{MAF},$gtinfo{$normalid}{DP},
		   sprintf("%.2f",abs($val{'tumorAltFrequency'} - $tmaf))),"\n";
  }
  foreach $c (keys %{$vars{$caseID}}) {
    foreach $p (keys %{$vars{$caseID}{$c}}) {
      next if $found{$c}{$p};
      %val = %{$vars{$caseID}{$c}{$p}};
      $val{'normalAltFrequency'} = 0 unless (exists $val{'normalAltFrequency'});
      print SNV join("\t",$caseID,$c,$p,$val{'geneName'},$val{'notation'},
		     sprintf("%.2f",$val{'tumorAltFrequency'}),
		     sprintf("%.2f",$val{'normalAltFrequency'})),"\n";
    }
  }
}
close SNV;

open OUT, ">cnv_validation.txt" or die $!;
print OUT join("\t",'CaseID','SampleID','Gene','Reported-CN',
	       'Reported-AbberationType','CHROM','Start','End',
	       'Reported-Cytoband','Score','CN','AbberationType',
	       'CHROM','Start','End','Cytoband','Score'),"\n";
foreach $caseID (keys %cnvs) {
  my @tumordirs = `ls -d $analdir/$caseID/*T_DNA_heme183`;
  chomp(@tumordirs);
  foreach $tdir (@tumordirs) {
    my ($dir,$pid,$tid) = split(/\//,$tdir);
    my $cnvfile = $dir."/".$pid."/".$tid."/".$tid.".cnv.answer.txt";
    unless (-e $cnvfile) {
	print $cnvfile,"\n";
	next;
    }
    open IN, "<$cnvfile" or die $!;
    my $head = <IN>;
    chomp($head);
    my @colnames = split(/\t/,$head);
    my %results;
    while (my $line = <IN>) {
      chomp($line);
      my @row = split(/\t/,$line);
      my %hash;
      foreach my $i (0..$#row) {
	$hash{$colnames[$i]} = $row[$i];
      }
      my $g = $hash{Gene};
      $results{$g} = 1;
      if ($cnvs{$caseID}{$g}) {
	print OUT join("\t",$caseID,$tid,$g,$cnvs{$caseID}{$g}{'copyNumber'},
		       $cnvs{$caseID}{$g}{'aberrationType'},
		       $cnvs{$caseID}{$g}{chrom},$cnvs{$caseID}{$g}{start},
		       $cnvs{$caseID}{$g}{end},$cnvs{$caseID}{$g}{cytoband},
		       $cnvs{$caseID}{$g}{score},$hash{CN},
		       $hash{'Abberation Type'},$hash{Chromosome},$hash{Start},
		       $hash{End},$hash{CytoBand},$hash{Score}),"\n";
      }
    }
    foreach $g (keys %{$cnvs{$caseID}}) {
      next if $results{$g};
      print OUT join("\t",$caseID,$tid,$g,$cnvs{$caseID}{$g}{'copyNumber'},
		     $cnvs{$caseID}{$g}{'aberrationType'},
		     $cnvs{$caseID}{$g}{chrom},$cnvs{$caseID}{$g}{start},
		     $cnvs{$caseID}{$g}{end},$cnvs{$caseID}{$g}{cytoband},
		     $cnvs{$caseID}{$g}{score}),"\n";
    }
  }
}
close OUT;
