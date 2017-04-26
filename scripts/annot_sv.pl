#!/usr/bin/perl -w
#svvcf2bed.pl

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my %opt = ();
my $results = GetOptions (\%opt,'input|i=s','refdb|r=s','help|h');

my $vcf = $opt{input};

my %lines;
my %eventid;
my $ct = 0;
open IN, "<$vcf" or die $!;
while (my $line = <IN>) {
  chomp($line);
  if ($line =~ m/^#CHROM/) {
    my @header = split(/\t/,$line);
    ($chrom, $pos,$id,$ref,$alt,$score,
     $filter,$info,$format,@subjacc) = split(/\t/, $line);
  }
  next if $line =~ m/#/;
  my ($chrom, $pos,$id,$ref,$alt,$score,
      $filter,$annot,$format,@gts) = split(/\t/, $line);
  $ct ++;
  my %hash = ();
  foreach $a (split(/;/,$annot)) {
    my ($key,$val) = split(/=/,$a);
    $hash{$key} = $val unless ($hash{$key});
  }
  if ($id eq 'N') {
      $id = 'NB'.sprintf("%06s",$ct);
  }
  my $evid = (split(/_/,$id))[0];
  $hash{'END'} = $pos+1 unless $hash{'END'};
  @deschead = split(":",$format);
 F1:foreach $sample (@subjacc) {
    my $allele_info = shift @gts;
    @ainfo = split(/:/, $allele_info);
    my %gtinfo = ();
    foreach $k (0..$#deschead) {
      $gtinfo{$deschead[$k]} = $ainfo[$k];
    }
    unless ($gtinfo{SU}) {
	$gtinfo{SU} = 0;
	$gtinfo{SU} = $gtinfo{RV}+$gtinfo{DV} if ($gtinfo{RV} && $gtinfo{DV});
    }
    $gtinfo{CN} = '' unless  $gtinfo{CN};
    if ($alt =~ m/chr(\w+):(\d+)/i) {
	if ($1 eq $chrom) {
	    my $locus = join(":",$chrom,$pos,$2);
	    $eventid{$id} = $locus;
	    $lines{$evid}{$locus}{$sample} =  join("\t",$hash{SVTYPE},$gtinfo{CN},$gtinfo{SU});
	}elsif ($id =~ m/_\d+/) {
	    my $locus = join(":",$chrom,$pos,$hash{END});
	    $eventid{$id} = $locus;
	    $lines{$evid}{$locus}{$sample} =  join("\t",$hash{SVTYPE},$gtinfo{CN},$gtinfo{SU});
	}else {
	    my $locus1 = join(":",$chrom,$pos,$hash{END});
	    my $id1 = $id."_1";
	    $eventid{$id1} = $locus1;
	    my $locus2 = join(":",'chr'.$1,$2,$2+1);
	    my $id2 = $id."_2";
	    $eventid{$id2} = $locus2;
	    $lines{$evid}{$locus1}{$sample} =  join("\t",$hash{SVTYPE},$gtinfo{CN},$gtinfo{SU});
	    $lines{$evid}{$locus2}{$sample} =  join("\t",$hash{SVTYPE},$gtinfo{CN},$gtinfo{SU});
	}
    }else {
	if ($hash{CHR2} && $hash{CHR2} eq $chrom) {
	    my $locus = join(":",$chrom,$pos,$hash{END});
	    $eventid{$id} = $locus;
	    $lines{$evid}{$locus}{$sample} =  join("\t",$hash{SVTYPE},$gtinfo{CN},$gtinfo{SU});
	}elsif ($hash{CHR2} && $hash{CHR2} ne $chrom) {
	    my $locus1 = join(":",$chrom,$pos,$pos+1);
	    my $id1 = $id."_1";
	    $eventid{$id1} = $locus1;
	    my $locus2 = join(":",$hash{CHR2},$hash{END},$hash{END}+1);
	    my $id2 = $id."_2";
	    $eventid{$id2} = $locus2;
	    $lines{$evid}{$locus1}{$sample} =  join("\t",$hash{SVTYPE},$gtinfo{CN},$gtinfo{SU});
	    $lines{$evid}{$locus2}{$sample} =  join("\t",$hash{SVTYPE},$gtinfo{CN},$gtinfo{SU});
	}unless ($hash{CHR2}) {
	    my $locus = join(":",$chrom,$pos,$hash{END});
	    $eventid{$id} = $locus;
	    $lines{$evid}{$locus}{$sample} =  join("\t",$hash{SVTYPE},$gtinfo{CN},$gtinfo{SU});
	}
    }
 }
}

close IN;

my @files = ("exonoverlap_sv.txt","geneoverlap_sv.txt");
foreach $file (@files) {
  open IN, "<$file" or die $!;
  while (my $line = <IN>) {
    chomp($line);
    next if ($line =~ m/#/);
    my ($chrom,$start,$end,$id,$chr,$start2,$end2,$gene) = split(/\t/, $line);
    my $evid = (split(/_/,$id))[0];
    my $locus = $eventid{$id};
    unless ($locus) {
	$locus = $eventid{$evid};
    }
    if ($file eq 'exonoverlap_sv.txt') {
	push @{$gene{$locus}}, $gene;
	$inc{$evid} = 1;
    }elsif ($file eq 'geneoverlap_sv.txt') {
	push @{$gene{$locus}}, $gene;
	$inc{$evid} = 1;
    }
  }
}

my $outfile = $vcf;
$outfile =~ s/vcf/annot.txt/g;

open OUT, ">$outfile" or die $!;
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
	  next unless $exonnum;
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
    foreach $sample (sort {$a cmp $b} keys %{$lines{$events}{$locus}}) {
      foreach $genesym (keys %gene_annots) {
	  my ($SVTYPE,$CN,$SU) = split(/\t/,$lines{$events}{$locus}{$sample});
	  next if $SU < 3;
	  print OUT join("\t",'SV.'.$events,split(":",$locus),$lines{$events}{$locus}{$sample},
			 $genesym, join(";",sort {$a cmp $b} keys %{$gene_annots{$genesym}}),$sample),"\n";
      }
    }
  }
}

