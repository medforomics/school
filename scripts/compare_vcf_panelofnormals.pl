#!/usr/bin/perl -w
#integrate_datasets.pl

#module load vcftools/0.1.14 samtools/1.6 bedtools/2.26.0 
open OM, "</project/shared/bicf_workflow_ref/GRCh38/clinseq_prj/panel1385.genelist.txt" or die $!;
while (my $line = <OM>) {
  chomp($line);
  $keep{$line} = 1;
}
close OM;
open OM, "</project/shared/bicf_workflow_ref/GRCh38/clinseq_prj/cancer.genelist.txt" or die $!;
while (my $line = <OM>) {
  chomp($line);
  $cgenelist{$line} = 1;
}
close OM;

my @vcflist = @ARGV;
foreach my $vfile (@vcflist) {
  open IN, "gunzip -c $vfile |" or die $!;
  my $tumorid;
 W1:while (my $line = <IN>) {
    chomp($line);
    if ($line =~ m/^#CHROM/) {
      my @header = split(/\t/,$line);
      ($chrom, $pos,$id,$ref,$alt,$score,
       $filter,$info,$format,@gtheader) = split(/\t/, $line);
      $tumorid = $gtheader[0];
      @sampids = @gtheader;
      next;
    }
    next if ($line =~ m/^#/);
    my ($chrom, $pos,$id,$ref,$alt,$score,
	$filter,$annot,$format,@gts) = split(/\t/, $line);
    next if ($ref =~ m/\./ || $alt =~ m/\./ || $alt=~ m/,X/);
    my %hash = ();
    foreach $a (split(/;/,$annot)) {
      my ($key,$val) = split(/=/,$a);
      $hash{$key} = $val unless ($hash{$key});
    }
    $hash{'HG38Loci'} = join(":",$chrom,$pos);
    my %fail;
    $fail{'UTSWBlacklist'} = 1 if ($hash{UTSWBlacklist});
    my @exacaf;
    if ($hash{AF_POPMAX}) {
	foreach (split(/,/,$hash{AF_POPMAX})) {
	    push @exacaf, $_ if ($_ ne '.');
	}
	@exacaf = sort {$b cmp $a} @exacaf;
    }
    $exacaf = $exacaf[0] if ($exacaf[0]);
    $fail{'COMMON'} = 1 unless ($exacaf && $exacaf <= 0.01);
    $fail{'StrandBias'} = 1 if (($hash{FS} && $hash{FS} > 60) || $filter =~ m/strandBias/i);
    my $cosmicsubj = 0;
    if ($hash{CNT}) {
	my @cosmicct = split(/,/,$hash{CNT}); 
	foreach $val (@cosmicct) {
	    $cosmicsubj += $val if ($val =~ m/^\d+$/);
      }
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
	$hash{$deschead[$k]} = $ainfo[$k] if ($subjid eq $tumorid);
      }
      $gtinfo{$subjid}{DP} = (split(/,/,$gtinfo{$subjid}{DP}))[0] if ($gtinfo{$subjid}{DP});
      next F1 unless ($gtinfo{$subjid}{DP} && $gtinfo{$subjid}{DP} ne '.' && $gtinfo{$subjid}{DP} >= 1);
      my @altct = split(/,/,$gtinfo{$subjid}{AD});
      my $refct = shift @altct;
      @altct2 = split(/,/,$gtinfo{$subjid}{AO});
      my $total = $refct;
      foreach  my $act (@altct) {
	next if ($act eq '.');
	$total += $act;
	push @mutallfreq, sprintf("%.4f",$act/$gtinfo{$subjid}{DP});
      }
    $gtinfo{$subjid}{MAF} = \@mutallfreq;
    }
    next unless ($gtinfo{$tumorid}{DP} && $gtinfo{$tumorid}{DP} ne '.' && $gtinfo{$tumorid}{DP} >= 20);
    unless ($gtinfo{$tumorid}{AO} =~ m/\d+/ && $gtinfo{$tumorid}{AD} =~ m/,/) {
      warn "Missing Alt:$line\n";
    }
    @tumormaf = @{$gtinfo{$tumorid}{MAF}};
    @tumoraltct = split(/,/,$gtinfo{$tumorid}{AO});
    
    if (exists $hash{INDEL}) {
      $hash{TYPE} = 'indel';
    }
    $hash{TYPE} = 'ambi' unless ($hash{"TYPE"});
    next if ($tumoraltct[0] eq '.');
    $hash{AF} = join(",",@tumormaf);
    my @callers ;
    @callers = split(/,/, $hash{CallSet}) if ($hash{CallSet});
    $hash{CallSet} = join(",",@callers);
    if ($id =~ m/COS/ && $cosmicsubj >= 5) {
      $fail{'LowAltCt'} = 1 if ($tumoraltct[0] < 3);
      $fail{'LowMAF'} = 1 if ($tumormaf[0] < 0.01);
    }else {
      $fail{'OneCaller'} = 1 if (scalar(@callers) < 2);
      $fail{'LowAltCt'} = 1 if ($tumoraltct[0] < 8);
      $fail{'LowMAF'} = 1 if ($tumormaf[0] < 0.05);
    }
    $freqtype = 'lowfreq';
    if ($tumormaf[0] > 0.30) {
	$type = 'germline';
    }
    my $newformat = 'GT:DP:AD:AO:RO';
    my @newgt;
    foreach $sample (@sampids) {
      my @gtdata;
      foreach $gt (split(/:/,$newformat)) {
	$gtinfo{$sample}{$gt} = '.' unless (exists $gtinfo{$sample}{$gt});
	push @gtdata, $gtinfo{$sample}{$gt};
      }
      push @newgt, join(":",@gtdata);
    }
    next unless ($hash{ANN});
    my $cancergene = 0;
    my $keepforvcf;
    foreach $trx (split(/,/,$hash{ANN})) {
      my ($allele,$effect,$impact,$gene,$geneid,$feature,
	$featureid,$biotype,$rank,$codon,$aa,$pos_dna,$len_cdna,
	$cds_pos,$cds_len,$aapos,$aalen,$distance,$err) = split(/\|/,$trx);
      next unless ($impact =~ m/HIGH|MODERATE/ || $effect =~ /splice/i);
      next unless $keep{$gene};
      $keepforvcf = $gene;
      $cancergene = 1 if ($cgenelist{$gene});
    }
    next unless $keepforvcf;
    my @fail = sort {$a cmp $b} keys %fail;
    if (scalar(@fail) == 0) {
      $filter = 'PASS';
    }else {
      $filter = join(";", 'FailedQC',@fail);
    }
    my @nannot;
    foreach $info (sort {$a cmp $b} keys %hash) {
      if (defined $hash{$info}) {
	push @nannot, $info."=".$hash{$info};
      }else {
	push @nannot, $info;
      }
    }
    $newannot = join(";",@nannot);
    $info{$chrom}{$pos}{$vfile} = [$freqtype,$hash{TYPE},$cancergene,$id,$ref,$alt,$newannot,$newformat,
				   @newgt] if ($filter eq 'PASS' && $id =~ /COSM/);
 }
  close IN;
}
open OUT, ">common_rare_snvs.txt" or die $!;
print OUT join("\t","Chr","Pos","NumberOfFiles","FileName","GermlineSomatic","SType","TierOneGene","ID",
	       "REF","ALT","Annot","Format","GTInfo"),"\n";
foreach $chrom (keys %info) {
    foreach $pos (keys %{$info{$chrom}}) {
	my @files = keys %{$info{$chrom}{$pos}};
	next if (scalar(@files) < 4);
	foreach $vf (@files) {
	    print OUT join("\t",$chrom,$pos,scalar(@files),$vf,@{$info{$chrom}{$pos}{$vf}}),"\n";
	}
    }
}
