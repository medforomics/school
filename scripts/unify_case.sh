#!/bin/bash
#union.sh

usage() {
  echo "-h Help documentation for gatkrunner.sh"
  echo "-r  --Reference Genome: GRCh38 or GRCm38"
  echo "-p  --Prefix for output file name"
  echo "-t  --TumorID"
  echo "-v  --tumor vcf"
  echo "-s  --somatic vcf"
  echo "-x  --rnaseq vcf"
  echo "-b  --rnaseq bam"
  echo "Example: bash union.sh -p prefix -r /path/GRCh38"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :r:p:t:n:v:s:x:b:f:h opt
do
    case $opt in
        r) index_path=$OPTARG;;
        p) subject=$OPTARG;;
        t) tumor_id=$OPTARG;;
        n) normal_id=$OPTARG;;
        v) tumor_vcf=$OPTARG;;
        s) somatic_vcf=$OPTARG;;
        x) rnaseq_vcf=$OPTARG;;
        b) rbam=$OPTARG;;
	f) fpkm=$OPTARG;;
        h) usage;;
    esac
done
if [[ -z $tumor_vcf ]] || [[ -z $subject ]] || [[ -z $index_path ]]; then
    usage
fi 

shift $(($OPTIND -1))
baseDir="`dirname \"$0\"`"
vepdir='/project/shared/bicf_workflow_ref/vcf2maf'

module load bedtools/2.26.0 samtools/1.6 vcftools/0.1.14
if [[ -a $somatic_vcf ]] 
then
    perl $baseDir/filter_somatic.pl $tumor_id $normal_id $somatic_vcf
    bgzip somatic.vcf
    tabix somatic.vcf.gz
    cp somatic.vcf.gz ${subject}.SS2.vcf.gz
    tabix -f $tumor_vcf
    bcftools annotate -Ov -a $tumor_vcf -o somatic.only.vcf --columns CHROM,POS,CallSet somatic.vcf.gz
    perl $baseDir/vcf2bed.pl somatic.only.vcf |cut -f 1,2,3 > somatic.bed
    bedtools intersect -header -v -b somatic.bed -a $tumor_vcf |perl -pe 's/\.final//g' > tumoronly.vcf
    vcf-shuffle-cols -t somatic.only.vcf tumoronly.vcf |bgzip > tumor.vcf.gz
    bgzip somatic.only.vcf
    vcf-concat somatic.only.vcf.gz tumor.vcf.gz |vcf-sort |bgzip > somatic_germline.vcf.gz
else
    ln -s $tumor_vcf somatic_germline.vcf.gz
fi

tabix somatic_germline.vcf.gz

if [[ -a $rbam ]]
then
    zcat somatic_germline.vcf.gz > dnavars.vcf
    perl $baseDir/vcf2bed.pl dnavars.vcf | perl -p -e 's/^chr//' |cut -f 1,2,3 > tumor.bed
    samtools index $rbam
    /project/shared/bicf_workflow_ref/seqprg/bam-readcount/bin/bam-readcount -w 0 -q 0 -b 25 -l tumor.bed -f ${index_path}/hisat_genome.fa $rbam > rnaseq.bamreadct.txt
fi

perl $baseDir/integrate_vcfs.pl ${subject} $tumor_id $normal_id $index_path $rnaseq_vcf
vcf-sort ${subject}.all.vcf | bedtools intersect -header -a stdin -b ${index_path}/UTSWV2.bed  |uniq |bgzip > ${subject}.utsw.vcf.gz
bgzip ${subject}.pass.vcf
tabix ${subject}.utsw.vcf.gz
tabix ${subject}.pass.vcf.gz
bcftools view -Oz -o ${subject}.vcf.gz -s ${tumor_id}  ${subject}.utsw.vcf.gz

#Makes TumorMutationBurenFile
bedtools intersect -header -a ${subject}.pass.vcf.gz -b ${index_path}/UTSWV2.bed  |uniq |bgzip > ${subject}.utswpass.vcf.gz
zgrep -c "SS=2" ${subject}.utswpass.vcf.gz |awk '{print "Class,TMB\n,"sprintf("%.2f",$1/4.6)}' > ${subject}.tmb.txt

#Convert to HG37
module load crossmap/0.2.5
CrossMap.py vcf ${index_path}/../hg38ToHg19.over.chain.gz ${subject}.utswpass.vcf.gz /project/apps_database/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa ${subject}.PASS.hg19.vcf
cat ${subject}.PASS.hg19.vcf |perl -p -e 's/^chr//g' > ${subject}.formaf.vcf

perl $baseDir/philips_excel.pl ${subject}.PASS.hg19.vcf $rnaseq_fpkm

#python $baseDir/../IntellispaceDemographics/gatherdemographics.py -i $subject -u phg_workflow -p $password -o ${subject}.xml
perl $baseDir/compareTumorNormal.pl ${subject}.utswpass.vcf.gz

zcat ${subject}.pass.vcf.gz |perl -p -e 's/^chr//g' > ${subject}.formaf.vcf
perl ${vepdir}/vcf2maf.pl --input ${subject}.formaf.vcf --output ${subject}.maf --species homo_sapiens --ncbi-build GRCh37 --ref-fasta ${vepdir}/.vep/homo_sapiens/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa --filter-vcf ${vepdir}/.vep/homo_sapiens/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz --cache-version 91 --vep-path ${vepdir}/variant_effect_predictor --tumor-id $tumor_id --normal-id $nromal_id --custom-enst ${index_path}/isoform_overrides_uniprot.txt --maf-center UTSW --vep-data ${vepdir}/.vep

#sequencestats
#exoncoverage
