#!/usr/bin/perl -w
#validation_json2txt.pl

use JSON;

my %keep;
open OM, "</project/shared/bicf_workflow_ref/human/grch38_cloud/panels/UTSW_V4_heme/genelist.txt" or die $!;
while (my $line = <OM>) {
    chomp($line);
    $keep{$line} = 1;
}

my @vtypes = ('variants','cnvs','translocations');
my %cnvs;
my %fusions;
my %vars;

foreach $jsonfile (@ARGV) {
  $jsontxt = `cat $jsonfile`;
  $jsonref  = decode_json($jsontxt);
  $prefix = (split(/\./,$jsonfile))[0];
  foreach my $vtype (@vtypes) {
    my %csam;
    foreach $tref (@{$jsonref->{$vtype}}) {
      my %hash = %{$tref};
      if ($vtype eq 'variants') {
	$vars{$prefix}{$hash{'geneName'}}{$hash{'notation'}} = \%hash;
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

open OUT, ">cnv_validation.txt" or die $!;
print OUT join("\t",'CaseID','Reported-CN','Reported-AbberationType','CHROM',
	       'Start','End','Reported-Cytoband','Score','CN','AbberationType',
	       'CHROM','Start','End','Cytoband','Score'),"\n";
foreach $caseID (keys %cnvs) {
  my @tumordirs = `ls -d $caseID/*T_DNA_pancancer1505`;
  chomp(@tumordirs);
  foreach $tdir (@tumordirs) {
    my ($pid,$tid) = split(/\//,$tdir);
    my $cnvfile = $pid."/".$tid."/".$tid.".cnv.answer.txt";
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

#print OUT join("\t",$prefix,$hash{chrom},$hash{pos},$hash{reference},
#	     $hash{alt},$hash{'tumorTotalDepth'},
#	     sprintf("%.2f",$hash{'tumorAltFrequency'}),
#	     $hash{'normalTotalDepth'},
#	     sprintf("%.2f",$hash{'normalAltFrequency'}),
#	     $hash{'geneName'},$hash{'notation'}),"\n";
#print OUT join("\t",$prefix,$hash{'caseId'},$gid,$hash{'copyNumber'},$hash{'aberrationType'},$hash{cytoband},$hash{score}),"\n";


