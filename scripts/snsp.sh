#!/bin/bash
#snsp.sh

usage() {
  echo "-h Help documentation for gatkrunner.sh"
  echo "-t  --Target Panel Bed File"
  echo "-p  --SampleID"
  echo "-b  --BAM File"
  echo "-v  --VCF File"
  echo "Example: bash snsp.sh -p prefix -b bamfile -t targetpanel -v VCF"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :t:r:v:s:p:b:h opt
do
    case $opt in
        p) id=$OPTARG;;
        t) tpanel=$OPTARG;;
	r) index_path=$OPTARG;;
        b) bam=$OPTARG;;
        v) vcf=$OPTARG;;
	h) usage;;
    esac
done

shift $(($OPTIND -1))
baseDir="`dirname \"$0\"`"

if [[ -z ${id} ]]; then
    usage
fi 
if [[ -z $index_path ]]; then
    index_path="/project/shared/bicf_workflow_ref/human/GRCh38/clinseq_prj"
fi
if [[ -z $index_path ]]; then
    tpanel="${index_path}/utswv2_cds.bed"
fi

shift $(($OPTIND -1))
baseDir="`dirname \"$0\"`"
source /etc/profile.d/modules.sh
module load samtools/gcc/1.8 htslib/gcc/1.8 bedtools/2.26.0 vcftools/0.1.14 bcftools/1.6

cut -f 1,2,3 ${index_path}/../gencode.CDS.bed |grep ^chr > gencode.cds.bed 

bedtools intersect -header -a ${index_path}/giab3.vcf.gz -b ${tpanel} | bedtools intersect -v -header -a stdin -b ${index_path}/HLA_HG38.bed | bedtools intersect -header -a stdin -b ${index_path}/giab3_platinum2.hg38.highConf.bed | bedtools intersect -header -a stdin -b gencode.cds.bed | vcf-sort | uniq | /project/shared/bicf_workflow_ref/seqprg/vt/vt decompose_blocksub - -a -o giab3.utswcoding.vcf
bgzip giab3.utswcoding.vcf
tabix giab3.utswcoding.vcf.gz
bedtools intersect -header -a ${index_path}/platinum_v2.vcf.gz -b ${tpanel} | bedtools intersect -v -header -a stdin -b ${index_path}/HLA_HG38.bed | bedtools intersect -header -a stdin -b ${index_path}/giab3_platinum2.hg38.highConf.bed  | bedtools intersect -header -a stdin -b gencode.cds.bed | vcf-sort | uniq | /project/shared/bicf_workflow_ref/seqprg/vt/vt decompose_blocksub - -a -o platinum_v2.utswcoding.vcf
bgzip platinum_v2.utswcoding.vcf
tabix platinum_v2.utswcoding.vcf.gz
bedtools intersect -header -a ${vcf} -b ${tpanel} | bedtools intersect -v -header -a stdin -b ${index_path}/HLA_HG38.bed | bedtools intersect -header -a stdin -b ${index_path}/giab3_platinum2.hg38.highConf.bed  | bedtools intersect -header -a stdin -b gencode.cds.bed | vcf-sort |uniq | bgzip > ${id}.utswcoding.vcf.gz
tabix ${id}.utswcoding.vcf.gz

bcftools annotate -Oz -a giab3.utswcoding.vcf.gz --columns CHROM,POS,REF,ALT,platforms -o ${id}.g.vcf.gz ${id}.utswcoding.vcf.gz
tabix ${id}.g.vcf.gz
bcftools annotate -Ov -a platinum_v2.utswcoding.vcf.gz --columns CHROM,POS,REF,ALT,MTD -o ${id}.eval.vcf ${id}.g.vcf.gz
bedtools genomecov -bga -split -ibam ${bam} -g ${index_path}/../genomefile.txt |awk '$4 > 19' | cut -f 1,2,3 | bedtools merge -i stdin > enoughcov.bed
bedtools intersect -header -a giab3.utswcoding.vcf.gz -b platinum_v2.utswcoding.vcf.gz |bedtools intersect -header -v -a stdin -b ${id}.utswcoding.vcf.gz | bedtools intersect -a stdin -b enoughcov.bed > fn.vcf
perl $baseDir/calc_snsp.pl ${id}.eval.vcf
