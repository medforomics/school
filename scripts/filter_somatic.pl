#!/usr/bin/perl -w
#use strict;

open OM, "</project/shared/bicf_workflow_ref/GRCh38/panel1385.genelist.txt" or die $!;
while (my $line = <OM>) {
    chomp($line);
    $keep{$line} = 1;
}
close OM;

my $prefix = shift @ARGV;
my $tid = shift @ARGV;
my $input = "$prefix\.annot.vcf.gz" or die $!;
open SOM, ">$prefix\.somatic.vcf" or die $!;
open IN, "gunzip -c $input|" or die $!;
 W1:while (my $line = <IN>) {
     chomp($line);
     if ($line =~ m/^#CHROM/) {
	 my @header = split(/\t/,$line);
	 ($chrom, $pos,$id,$ref,$alt,$score,
	  $filter,$info,$format,@gtheader) = split(/\t/, $line);
     }
     if ($line =~ m/^#/) {
	 print  SOM $line,"\n";
	 next;
     }
     my ($chrom, $pos,$id,$ref,$alt,$score,
	 $filter,$annot,$format,@gts) = split(/\t/, $line);
     next if ($ref =~ m/\./ || $alt =~ m/\./ || $alt=~ m/,X/);
     my %hash = ();
     foreach $a (split(/;/,$annot)) {
	 my ($key,$val) = split(/=/,$a);
	 $hash{$key} = $val unless ($hash{$key});
     }
     my %fail;
     $fail{'UTSWBlacklist'} = 1 if ($hash{UTSWBlacklist});
     my $exacaf = '';
     if ($hash{AC_POPMAX} && $hash{AN_POPMAX}) {
	 @exacs = split(/,/,$hash{AC_POPMAX});
	 my $ac = 0;
	 foreach $val (@exacs) {
	     $ac += $val if ($val =~ m/^\d+$/);
	 }
	 @exans = split(/,/,$hash{AN_POPMAX});
	 my $an = 0;
	 foreach $val (@exans) {
	     $an += $val if ($val =~ m/^\d+$/);
	 }
	 $exacaf = sprintf("%.4f",$ac/$an) if ($ac > 0 && $an > 10);
     }
     unless ($exacaf eq '' || $exacaf <= 0.01) {
	 #$fail{'COMMON'} = 1;
     }
     my $cosmicsubj = 0;
     if ($hash{CNT}) {
	 my @cosmicct = split(/,/,$hash{CNT}); 
	 foreach $val (@cosmicct) {
	     $cosmicsubj += $val if ($val =~ m/^\d+$/);
	 }
     }
     my @maf;
     my @dp;
     my @ao;
     my @genotypes = @gts;
     my @deschead = split(/:/,$format);
   F1:foreach my $subjid (@gtheader) {
       my $allele_info = shift @gts;
       @ainfo = split(/:/, $allele_info);
       my %gtinfo = ();
       my @mutallfreq = ();
       foreach my $k (0..$#ainfo) {
	   $gtinfo{$deschead[$k]} = $ainfo[$k];
       }
       next W1 if ($gtinfo{DP} < 10);
       my @altct = split(/,/,$gtinfo{AO});
       foreach  my $act (@altct) {
	   push @mutallfreq, sprintf("%.4f",$act/$gtinfo{DP});
       }
       push @dp, $gtinfo{DP};
       push @maf, \@mutallfreq;
       my @sortao = sort {$b <=> $a} @altct;
       push @ao, $sortao[0];
   }
     if ($gtheader[1] eq $tid) {
	 @maf = reverse(@maf);
	 @dp = reverse(@dp);
	 @ao = reverse(@ao);
	 @genotypes = reverse(@genotypes);
     }
     $hash{AF} = join(",",@{$maf[0]});
     $hash{NormalAF} =  join(",",@{$maf[1]});
     $hash{DP} = $dp[0];
     $hash{NormalDP} = $dp[1];
     next if ($maf[1][0] > 0.005 || $maf[1][0]*5 > $maf[0][0]);
     my $newgt = $genotypes[0];
     foreach (@dp) {
	 $fail{'LowDepth'} = 1 if ($_ < 20);
    }
     my @callers = split(/,/,$hash{CallSet});
     if ($id =~ m/COS/ && $cosmicsubj >= 5) {
	 $fail{'LowAltCt'} = 1 if ($ao[0] < 3);
	 $fail{'LowMAF'} = 1 if ($maf[0][0] < 0.01);
     }else {
	 $fail{'OneCaller'} = 1 if (scalar(@callers) < 2);
	 $fail{'LowAltCt'} = 1 if ($ao[0] < 8);
	 $fail{'LowMAF'} = 1 if ($maf[0][0] < 0.05);
     }
     if ($rnaval{$chrom}{$pos}) {
	 $hash{RnaSeqValidation} = 1;
     } 
     $hash{Somatic} = 1;
     $hash{SomaticCallSet}=$hash{CallSet};
     next unless ($hash{ANN});
     foreach $trx (split(/,/,$hash{ANN})) {
	 my ($allele,$effect,$impact,$gene,$geneid,$feature,
	     $featureid,$biotype,$rank,$codon,$aa,$pos_dna,$len_cdna,
	     $cds_pos,$cds_len,$aapos,$aalen,$distance,$err) = split(/\|/,$trx);
	 next unless $keep{$gene};
	 if ($rnaseqct{$gene} && $rnaseqct{$gene} > 10) {
	     $hash{logcpm}=sprintf("%.1f",log2(1000000*$rnaseqct{$gene}/$rnaseqct{'total'}));
	 } if ($fpkm{$gene}) {
	     $hash{fpkm} = sprintf("%.1f",$fpkm{$gene});
	 }
	 my @fail = keys %fail;
	 if (scalar(@fail) == 0) {
	     $filter = 'PASS';
	     print SOM join("\t",$chrom, $pos,$id,$ref,$alt,$score,
			    $filter,$annot,$format,$newgt),"\n";
	 }else {
	     next W1;
	 }
     }
}
close IN;
