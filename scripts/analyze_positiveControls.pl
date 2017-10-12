#!/usr/bin/perl -w
use strict; 

#perl analyze_positiveControls --annot GM12878.annot.vcf.gz

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %opt = ();
my $results = GetOptions (\%opt,'annot|a=s','help|h');

my $prefix = $opt{annot};
$prefix =~ s/\.annot\.vcf\.gz//;
system("perl /project/PHG/PHG_Clinical/clinseq_workflows/scripts/filter_germline.pl $prefix");

system("bedtools intersect -header -a $prefix\.germline.vcf -b /project/shared/bicf_workflow_ref/GRCh38/utswv2_cds.bed | bedtools intersect -v -header -a stdin -b /project/shared/bicf_workflow_ref/GRCh38/HLA_HG38.bed | bedtools intersect -header -a stdin -b /project/shared/bicf_workflow_ref/GRCh38/giab_v2_highConf.bed | bgzip > $prefix\.utswcoding.vcf.gz");

system("bedtools multiinter -i $prefix\.utswcoding.vcf.gz /project/shared/bicf_workflow_ref/GRCh38/giab_v2.utsw.highConf.noHLA.vcf.gz /project/shared/bicf_workflow_ref/GRCh38/platinum_v2.utsw.highConf.noHLA.vcf.gz -names union giab platinum |cut -f 1,2,3,5 | bedtools sort -i stdin | bedtools merge -c 4 -o distinct | bgzip > $prefix\.utswcoding.multiinter.bed.gz");

system("tabix $prefix\.utswcoding.multiinter.bed.gz");

system("bcftools annotate -a  $prefix\.utswcoding.multiinter.bed.gz --columns CHROM,FROM,TO,PlatRef -h /project/shared/bicf_workflow_ref/GRCh38/PlatRef.header $prefix\.utswcoding.vcf.gz > $prefix\.eval.vcf");

system("perl /project/PHG/PHG_Clinical/clinseq_workflows/scripts/calc_snps.pl $prefix\.eval.vcf >$prefix\.calcsnps.txt");
