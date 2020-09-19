#!/usr/bin/perl -w
#extract_greyzone.pl

open GLIST, "</project/shared/bicf_workflow_ref/human/grch38_cloud/panels/UTSW_V4_pancancer/genelist.txt";
while (my $line = <GLIST>) {
    chomp($line);
    $keep{$line} = 1;
}

open REF, "</project/shared/bicf_workflow_ref/human/grch38_cloud/giab_ref/HD827_variant_annotation_COSMICv85_dbSNPv151.txt" or die $!;
$header = <REF>;
chomp($header);
my @colnames = split(/\t/,$header);
while (my $line = <REF>) {
    chomp($line);
    my @row = split(/\t/,$line);
    my %hash;
    foreach my $i (0..$#row) {
	$hash{$colnames[$i]} = $row[$i];
    }
    next unless ($keep{$hash{Gene}});
    $hash{Notes} = 'NA' unless ($hash{Notes});
    $missing{$hash{CHROM}}{$hash{POS}} = \%hash;
    $info{$hash{CHROM}}{$hash{POS}} = [$hash{Ref},$hash{Alt},$hash{'T-OAF'},$hash{'Notes'}];
}

open OUT, ">HD827_sn.csv" or die $!;
print OUT join(",",'Sample','LociGRCh38','ID',
	       'Gene','Nucleotide','AminoAcid','Effect','Ref','Alt',
	       'T-EAF','T-OAF','T-DP','Notes'),"\n";

my @vcffiles = @ARGV;
foreach $vcf (@vcffiles) {
  open VCF, "gunzip -c $vcf |" or die $!;
  $vcf =~ m/(HD\d+)/;
  $sampid = $1;
  while (my $line = <VCF>) {
    chomp($line);
    if ($line =~ m/^#CHROM/) {
      my @header = split(/\t/,$line);
      ($chrom, $pos,$id,$ref,$alt,$score,
       $filter,$info,$format,@gtheader) = split(/\t/, $line);
    }
    if ($line =~ m/^#/) {
      next;
    }
    my ($chrom, $pos,$id,$ref,$alt,$score,
	$filter,$annot,$format,@gts) = split(/\t/, $line);
    next unless  $info{$chrom}{$pos};
    my ($eref,$ealt,$teaf,$notes) = @{$info{$chrom}{$pos}};
    next unless ($ref eq $eref && $alt eq $ealt);
    next if ($ref =~ m/\./ || $alt =~ m/\./ || $alt=~ m/,X/);
    my %hash = ();
    foreach $a (split(/;/,$annot)) {
      my ($key,$val) = split(/=/,$a);
      #$val =~ s/,/\|/g if ($val);
      $hash{$key} = $val unless ($hash{$key});
    }
    my @deschead = split(/:/,$format);
    my $allele_info = shift @gts;
    my @ainfo = split(/:/, $allele_info);
    my @mutallfreq = ();
    foreach my $k (0..$#ainfo) {
      $gtinfo{$deschead[$k]} = $ainfo[$k];
    }
    $hash{DP} = (split(/,/,$gtinfo{DP}))[0];
    next F1 unless ($hash{DP} && $hash{DP} ne '.' && $hash{DP} >= 1);
    @altct = split(/,/,$gtinfo{AO});
    foreach  my $act (@altct) {
      next if ($act eq '.');
      push @mutallfreq, sprintf("%.4f",$act/$hash{DP});
    }
    $hash{AF} = \@mutallfreq;
    my $aachange;
    my $genechange;
    my $effectchange;
    my $aa;
    foreach $trx (split(/,/,$hash{ANN})) {
      my ($allele,$effect,$impact,$gene,$geneid,$feature,
	  $featureid,$biotype,$rank,$codon,$aa,$pos_dna,$len_cdna,
	  $cds_pos,$cds_len,$aapos,$aalen,$distance,$err) = split(/\|/,$trx);
      $genechange = $gene unless $genechange;
      $effectchange = $effect unless $effectchange;
      unless ($aachange) {
	$aachange = $trx if ($aa);
      }
    }
    if ($aachange) {
      ($allele,$effect,$impact,$gene,$geneid,$feature,
       $featureid,$biotype,$rank,$codon,$aa,$pos_dna,$len_cdna,
       $cds_pos,$cds_len,$aapos,$aalen,$distance,$err) = split(/\|/,$aachange);
    }else {
      $gene = $genechange;
      $effect = $effectchange;
    }
    $aa = 'NA' unless $aa;
    print OUT join(",",$sampid,join(":",$chrom,$pos),$id,
		   $gene,$codon,$aa,$effect,$ref,$alt,$teaf,
		   join(";",@{$hash{AF}}),$hash{DP},$notes),"\n";
    delete $missing{$chrom}{$pos};
  }
}

foreach $chrom (keys %missing) {
    foreach $pos (keys %{$missing{$chrom}}) {
	print OUT join("\t",$chrom,$pos,$missing{$chrom}{$pos}{Gene},
		   $missing{$chrom}{$pos}{COSM90AA},
		   $missing{$chrom}{$pos}{'T-OAF'},
		   $missing{$chrom}{$pos}{Notes}),"\n";
		   
    }
}
