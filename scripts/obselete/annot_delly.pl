#!/usr/bin/perl -w
#svvcf2bed.pl

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my %opt = ();
my $results = GetOptions (\%opt,'input|i=s','refdb|r=s','help|h');

my $vcf = $opt{input};

my %lines;
my %eventid;
my $sample;
open BED, ">$vcf\.bed" or die $!;
open IN, "gunzip -c $vcf|" or die $!;
while (my $line = <IN>) {
  chomp($line);
  if ($line =~ m/^#CHROM/) {
    my @header = split(/\t/,$line);
    ($chrom, $pos,$evid,$ref,$alt,$score,
     $filter,$info,$format,@subjacc) = split(/\t/, $line);
     $sample = $subjacc[0];

  }
  next if $line =~ m/#/;
  my ($chrom, $pos,$evid,$ref,$alt,$score,
      $filter,$annot,$format,@gts) = split(/\t/, $line);
  #if ($id eq 'TRA00000006') {
  #    warn "debugging\n";
  #}
  my %hash = ();
  foreach $a (split(/;/,$annot)) {
    my ($key,$val) = split(/=/,$a);
    $hash{$key} = $val unless ($hash{$key});
  }
  @deschead = split(":",$format);
  my $allele_info = shift @gts;
  @ainfo = split(/:/, $allele_info);
  my %gtinfo = ();
  foreach $k (0..$#deschead) {
      $gtinfo{$deschead[$k]} = $ainfo[$k];
  }
  unless ($gtinfo{SU}) {
      $gtinfo{SU} = 0;
      $gtinfo{SU} = $gtinfo{RV}+$gtinfo{DV} if (exists $gtinfo{RV} && exists $gtinfo{DV});
  }
  $gtinfo{CN} = '' unless  $gtinfo{CN};
  next if ($gtinfo{SU} < 2);
  $id = $evid;
  if ($hash{CHR2} ne $chrom) {
      print BED join("\t",$hash{CHR2},$hash{'END'},$hash{'END'}+1,$id."_2"),"\n";
      print BED join("\t",$chrom,$pos,$pos+1,$id."_2"),"\n";
      $locus1 = join(":",$chrom,$hash{'END'},$hash{'END'}+1);
      $locus2 = join(":",$chrom,$pos,$pos+1);
      $eventid{$id."_1"} = $locus2;
      $eventid{$id."_2"} = $locus1;
      $lines{$evid}{$locus1} =  join("\t",$hash{SVTYPE},$gtinfo{CN},$gtinfo{SU});
      $lines{$evid}{$locus2} =  join("\t",$hash{SVTYPE},$gtinfo{CN},$gtinfo{SU});
  }else {
      if ($hash{'END'}) {
	  if ($pos > $hash{'END'}) {
	      $temp = $pos;
	      $pos = $hash{'END'};
	      $hash{'END'} = $temp;
	  }
	  next if ($hash{'END'}-$pos < 500);
      }else {
	  $hash{'END'} = $pos+1;
      }
      print BED join("\t",$chrom,$pos,$hash{'END'},$id),"\n";
      $locus = join(":",$chrom,$pos,$hash{'END'});
      $eventid{$id} = $locus;
      $lines{$evid}{$locus}{$sample} =  join("\t",$hash{SVTYPE},$gtinfo{CN},$gtinfo{SU});
  }
}
close BED;
close IN;

system(qq{bedtools intersect -header -wb -v -a $vcf\.bed -b $opt{refdb}\/rmsk.bed > rmskoverlap_sv.txt});
system(qq{bedtools intersect -header -wb -a rmskoverlap_sv.txt -b $opt{refdb}\/dgv.bed > dgvoverlap_sv.txt});
system(qq{bedtools intersect -header -wb -a rmskoverlap_sv.txt -b $opt{refdb}\/gencode.exons.bed > exonoverlap_sv.txt});
system(qq{bedtools intersect -v -header -wb -a $vcf\.bed -b $opt{refdb}\/gencode.exons.bed | bedtools intersect -header -wb -a stdin -b $opt{refdb}\/gencode.genes.chr.bed > geneoverlap_sv.txt});

my @files = ("dgvoverlap_sv.txt","exonoverlap_sv.txt","geneoverlap_sv.txt");
foreach $file (@files) {
  open IN, "<$file" or die $!;
  while (my $line = <IN>) {
    chomp($line);
    next if ($line =~ m/#/);
    my ($chrom,$start,$end,$id,$chr,$start2,$end2,$gene) = split(/\t/, $line);
    my $evid = (split(/_/,$id))[0];
    my $locus = $eventid{$id};
    if ($file eq 'exonoverlap_sv.txt') {
	push @{$gene{$locus}}, $gene;
	$inc{$evid} = 1;
    }elsif ($file eq 'geneoverlap_sv.txt') {
	push @{$gene{$locus}}, $gene;
	$inc{$evid} = 1;
    }else {
	my $change = join("|",(split(/\|/,$gene))[0]);
	$dgv{$locus}{$change} = 1;
    }
  }
}
open OUT, ">$vcf\.annot.txt" or die $!;
foreach $events (sort keys %lines) {
    next unless ($inc{$events});
    foreach $locus (sort {$a cmp $b} keys %{$lines{$events}}) {
	my %gene_annots;
	if ($gene{$locus}) {
	    my @genes = @{$gene{$locus}};
	    my %hash;
	    foreach $g (@genes) {
		foreach $trx (split(/,/,$g)) {
		    my ($symbol,$trxid,$exonnum) = split(/\|/,$trx);
		    $exonnum = '' unless $exonnum;
		    $hash{$symbol}{$trxid}{$exonnum} = 1;
		} 
	    }
	    foreach $sym (keys %hash) {
		foreach $trxid (keys %{$hash{$sym}}) {
		    my @exons = sort {$a <=> $b} keys %{$hash{$sym}{$trxid}};
		    my $exon_bound = $exons[0];
		    if ($#exons > 0) {
			$exon_bound = join("-",$exons[0],$exons[-1]);
		    }
		    $gene_annots{$sym}{$exon_bound} = 1;
		}
	    }
	}
	my $dgv = '';
	$dgv = join(";",keys %{$dgv{$locus}}) if ($dgv{$locus});
	foreach $genesym (keys %gene_annots) {
	    print OUT join("\t",$events,split(":",$locus),$lines{$events}{$locus}{$sample},
			   $genesym, join(";",keys %{$gene_annots{$genesym}}),$sample,$dgv),"\n";
	}
    }
}


