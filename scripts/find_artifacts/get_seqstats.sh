#!/bin/bash

module load repeatModeler gatk/4.1.2.0 picard/2.10.3 htslib/gcc/1.8 samtools/gcc/1.8 bedtools/2.26.0 bcftools/gcc/1.8 snpeff/4.3q 
caseID=$1
baseDir="`dirname \"$0\"`"
cd /archive/PHG/PHG_Clinical/toarchive/caseDirs
if [[ ! -f ${caseID}.json ]]
then
    curl -o ${caseID}.json "https://nuclia.biohpc.swmed.edu/getSampleAndRunIds?token=$nucliatoken&projectId=$caseID"
fi

genoDir='/project/shared/bicf_workflow_ref/human/GRCh38/'
jsonout=$(perl ${baseDir}/parse_getsampleJson.pl ${caseID}.json 2>&1)
myarray=($jsonout)
tumorID=${myarray[3]}
dnarunid=${myarray[2]}

bam="${caseID}/${tumorID}/${tumorID}.consensus.bam"

samtools index -@ 4 $bam
zgrep -v "#" ${caseID}/${caseID}.vcf.gz | awk '{if(length($4) > length($5)) print $1"\t"($2-1)"\t"($2+length($4)-1)"\t"$1":"$2":"$4":"$5; else print $1"\t"($2-1)"\t"($2+length($5)-1)"\t"$1":"$2":"$4":"$5}' > ${caseID}/${caseID}.bed
awk '{print $1"\t"$2-25"\t"$3+25"\t"$4}' ${caseID}/${caseID}.bed > ${caseID}/${caseID}.50window.bed
awk '{print $1"\t"$2-5"\t"$3+5"\t"$4}' ${caseID}/${caseID}.bed > ${caseID}/${caseID}.10window.bed
bedtools intersect -a ${caseID}/${caseID}.10window.bed -b ${genoDir}/repeat_regions.bed.gz -wao > ${caseID}/${caseID}.repeat.vcf.txt
bedtools getfasta -fi ${genoDir}/genome.fa -bed ${caseID}/${caseID}.50window.bed -fo ${caseID}/${caseID}.50window.fasta
/project/shared/bicf_workflow_ref/seqprg/bin/seqkit fx2tab --gc ${caseID}/${caseID}.50window.fasta > ${caseID}/${caseID}.50window.gc.txt
dustmasker -in ${caseID}/${caseID}.50window.fasta -out ${caseID}/${caseID}.50window.dust.txt
/project/shared/bicf_workflow_ref/seqprg/bin/seqkit locate -f motifs.fa -o ${caseID}/${caseID}.50window.homopolymer.txt ${caseID}/${caseID}.50window.fasta
/project/shared/bicf_workflow_ref/seqprg/bamutils/bamUtil-1.0.14/bin/bam stats --regionList ${caseID}/${caseID}.bed --in $bam --pBaseQC ${caseID}/${caseID}.bamstat
/project/shared/bicf_workflow_ref/seqprg/bam-readcount/bin/bam-readcount -l ${caseID}/${caseID}.10window.bed -w 0 -q 0 -f /project/shared/bicf_workflow_ref/human/GRCh38/genome.fa $bam > ${caseID}/${caseID}.bamreadct.txt

list1=`ls ${caseID}/dna*/*ori.vcf.gz ${caseID}/dna*/*strelka2.vcf.gz ${caseID}/somatic*/*ori.vcf.gz |grep -v strelka2.ori | grep -v platypus`

for i in $list1; do
    prefix="${i%.vcf.gz}"
    tabix -f $i
    gatk VariantAnnotator -R /project/shared/bicf_workflow_ref/human/GRCh38/genome.fa -V $i -I $bam -A BaseQuality -A FisherStrand -A FragmentLength -A MappingQuality -A RMSMappingQuality -A ReadPosition -A BaseQualityRankSumTest -A ClippingRankSumTest -O ${prefix}.stat.vcf.gz
    tabix -f ${prefix}.stat.vcf.gz
    bcftools norm -m - ${prefix}.stat.vcf.gz | bedtools intersect -header -b ${caseID}/${caseID}.bed -a stdin | java -jar $SNPEFF_HOME/SnpSift.jar extractFields - CHROM POS REF ALT FS MBQ MFRL MMQ MPOS MQ BaseQRankSum ClippingRankSum > ${prefix}.stat.txt
done

perl ${baseDir}/merge_quals.pl $caseID $tumorID $dnarunid
