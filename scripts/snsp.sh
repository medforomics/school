#!/bin/bash
#snsp.sh

ID=$1
Bam=$2
TPANEL=$3

usage() {
  echo "-h Help documentation for gatkrunner.sh"
  echo "-t  --Target Panel Bed File"
  echo "-p  --SampleID"
  echo "-r  --Reference Genome: GRCh38 or GRCm38"
  echo "Example: bash snsp.sh -p prefix -b bamfile -t targetpanel"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :r:p:t:b:h opt
do
    case $opt in
        r) index_path=$OPTARG;;
        p) ID=$OPTARG;;
        t) targetbed=$OPTARG;;
        h) usage;;
    esac
done
if [[ -z $ID ]]; then
    usage
fi 
if [[ -z $index_path ]]; then
    index_path="/project/shared/bicf_workflow_ref/GRCh38/clinseq_prj"
fi
if [[ -z $index_path ]]; then
    targetbed="${index_path}/utswv2_cds.bed"
fi

shift $(($OPTIND -1))
baseDir="`dirname \"$0\"`"

module load samtools/1.6 bedtools/2.26.0 vcftools/0.1.14

perl $baseDir/filter_giab.pl $ID

#cat ${index_path}/../gencode.CDS.bed |awk '{print $1"\t"$2-5"\t"$3+5}' |grep ^chr > gencode.cds.bed 
cut -f 1,2,3 ${index_path}/../gencode.CDS.bed |grep ^chr > gencode.cds.bed 

bedtools intersect -header -a ${index_path}/giab3.vcf.gz -b ${targetbed} | bedtools intersect -v -header -a stdin -b ${index_path}/HLA_HG38.bed | bedtools intersect -header -a stdin -b ${index_path}/giab3_platinum2.hg38.highConf.bed  | bedtools intersect -header -a stdin -b gencode.cds.bed | vcf-sort |uniq |bgzip > giab3.utswcoding.vcf.gz
bedtools intersect -header -a ${index_path}/platinum_v2.vcf.gz -b ${targetbed} | bedtools intersect -v -header -a stdin -b ${index_path}/HLA_HG38.bed | bedtools intersect -header -a stdin -b ${index_path}/giab3_platinum2.hg38.highConf.bed  | bedtools intersect -header -a stdin -b gencode.cds.bed | vcf-sort |uniq | bgzip > platinum_v2.utswcoding.vcf.gz
bedtools intersect -header -a ${ID}.PASS.vcf -b ${targetbed} | bedtools intersect -v -header -a stdin -b ${index_path}/HLA_HG38.bed | bedtools intersect -header -a stdin -b ${index_path}/giab3_platinum2.hg38.highConf.bed  | bedtools intersect -header -a stdin -b gencode.cds.bed | vcf-sort |uniq | bgzip > ${ID}.utswcoding.vcf.gz
tabix ${ID}.utswcoding.vcf.gz

bedtools multiinter -i ${ID}.utswcoding.vcf.gz giab3.utswcoding.vcf.gz platinum_v2.utswcoding.vcf.gz -names union giab platinum |cut -f 1,2,3,5 | bedtools sort -i stdin | bedtools merge -c 4 -o distinct | bgzip > ${ID}.utswcoding.multiinter.bed.gz 
tabix ${ID}.utswcoding.multiinter.bed.gz
bcftools annotate -a  ${ID}.utswcoding.multiinter.bed.gz --columns CHROM,FROM,TO,PlatRef -h /project/shared/bicf_workflow_ref/PlatRef.header ${ID}.utswcoding.vcf.gz > ${ID}.eval.vcf
perl $baseDir/calc_snsp.pl ${ID}.eval.vcf
