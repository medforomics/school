#!/usr/bin/perl -w
#integrate_datasets.pl
#module load vcftools/0.1.14 samtools/1.6 bedtools/2.26.0 
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my @keep = ('200000001','200000007','338927776','340042036','340170358','340170802','340281254','340290745','340301098','340545456','341037908','341038022','341038233','341234023','341642139','341739113','341892701','342338180','342839359','343452345','345077679','345393749','346266256','346332321','346489686','347270798','347395456','347497830','347803419','347824061','348794361','349482726','350400307','350890705','352051477','352550973','352602958','352770619','354735429','356887943','358489369','360764309','360766818','360888243','361229076','361941978','361969893','362944875','363311696','363453532','364096698','364319588','364350014','365138718','365328489','365339853','70A0B0011','71Z83NF3P','730FHG980','754Q9J5M4');

my %keep = map {$_=>1} @keep;

my %bl;
open BL, "<blacklist.txt" or die $!;
while (my $line = <BL>) {
  chomp($line);
  $key = join(":",split(/\t/,$line));
  $bl{$key} = 1;
}
my %sampkeep;
open NOR, "<hasnormal.txt" or die $!;
while (my $line = <NOR>) {
  chomp($line);
  next unless $keep{$line};
  $sampkeep{$line} = 1;
}
my %opt = ();

print join("\t",'CaseID','WithNormal','WithoutNormal','Germline','Somatic','RNAseq'),"\n";

my @vcffiles = @ARGV;
foreach $vcf (@vcffiles) {
  open IN, "gunzip -c $vcf |" or die $!;
  my ($caseid,@ext) = split(/\./,$vcf);
  my ($ordid,$mrn) = split(/-/,$caseid);
  next unless ($sampkeep{$ordid});
  my $somatic_tmb = 0;
  my $tonly_tmb = 0;
  my $total = 0;
  my $germline = 0;
  my $rnaseq = 0;
 W1:while (my $line = <IN>) {
    chomp($line);
    if ($line =~ m/^#/) {
      next;
    }
    my ($chrom, $pos,$id,$ref,$alt,$score,
	$filter,$annot,$format,@gts) = split(/\t/, $line);
    next if ($ref =~ m/\./ || $alt =~ m/\./ || $alt=~ m/,X/);
    $key = join(":",$chrom,$pos,$ref,$alt);
    my %hash = ();
    foreach $a (split(/;/,$annot)) {
      my ($key,$val) = split(/=/,$a);
      $hash{$key} = $val unless ($hash{$key});
    }
    next unless ($filter eq 'PASS');
    $total ++;
    if ($hash{SS} == 1) {
	$germline ++;
    }
    if ($hash{SS} == 2) {
      $somatic_tmb ++;
    }
    if ($hash{RnaSeqValidation}) {
	$rnaseq ++;
    }
    next if ($hash{TYPE} eq 'indel');
    next if ($hash{ANN} =~ m/MODERATE|HIGH/);
    next if ($hash{AF_POPMAX} && $hash{AF_POPMAX} ne '.' && $hash{AF_POPMAX} > 0.0001);
    next if ($ref eq 'C' && $alt eq 'T' && $hash{AF} < 0.2);
    next if ($ref eq 'G' && $alt eq 'A' && $hash{AF} < 0.2);
    next if ($bl{$key});
    $tonly_tmb ++;
  }
  print join("\t",$ordid,$somatic_tmb,$tonly_tmb,
	     sprintf("%.4f",$germline/$total),
	     sprintf("%.4f",$somatic_tmb/$total),
	     sprintf("%.4f",$rnaseq/$total)),"\n";
}
close IN;
