#!/usr/bin/perl -w
#extract_greyzone.pl

open REF, "</project/shared/bicf_workflow_ref/human/GRCh38/clinseq_prj/HorizonReference.txt" or die $!;
while (my $line = <REF>) {
    chomp($line);
    my ($sample,$cosmid,$chr,$pos) = split(/\t/,$line);
    $loci = join(":",$chr,$pos);
    $info{$sample}{$loci} = $cosmid;
}

open OUT, ">Horizon.csv" or die $!;
print OUT join(",",'Sample','LociGRCh38','ID',
	       'Gene','Nucleotide','AminoAcid','Effect','Ref','Alt',
	       'Tumor DNA AF','Tumor DNA Depth','CallSet'),"\n";

my @vcffiles = @ARGV;
foreach $vcf (@vcffiles) {
  open VCF, "gunzip -c $vcf |" or die $!;
  $vcf =~ m/_(HD\d+)_/;
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
    $loci = join(":",$chrom,$pos);
    next unless  $info{$sampid}{$loci};
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
    print OUT join(",",$sampid,join(":",$chrom,$pos),$id,
		   $gene,$codon,$aa,$effect,$ref,$alt,
		   join(";",@{$hash{AF}}),$hash{DP},$hash{CallSet}),"\n";
  }
}

