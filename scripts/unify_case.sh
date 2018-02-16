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
while getopts :r:p:t:v:s:x:n:b:h opt
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
    echo $normal $tumor $algo
    usage
fi 

shift $(($OPTIND -1))
baseDir="`dirname \"$0\"`"

module load bedtools/2.26.0 samtools/1.6 vcftools/0.1.14
if [[ -a $somatic_vcf ]] 
then
    perl $baseDir/filter_somatic.pl $tumor_id $normal_id $somatic_vcf
    bgzip somatic.vcf
    tabix somatic.vcf.gz
    tabix $tumor_vcf
    bcftools annotate -Ov -a $tumor_vcf -o somatic.only.vcf --columns CHROM,POS,CallSet somatic.vcf.gz
    perl $baseDir/vcf2bed.pl somatic.only.vcf |cut -f 1,2,3 > somatic.bed
    bedtools intersect -header -v -b somatic.bed -a $tumor_vcf |perl -p -e 's/\.final//g' > tumoronly.vcf
    vcf-shuffle-cols -t somatic.only.vcf tumoronly.vcf |bgzip > tumor.vcf.gz
    bgzip somatic.only.vcf
    vcf-concat somatic.only.vcf.gz tumor.vcf.gz |vcf-sort |bgzip > somatic_germline.vcf.gz
else
    ln -s $tumor_vcf somatic_germline.vcf.gz
fi
tabix somatic_germline.vcf.gz
if [[ -a $rnaseq_vcf ]]
then
    zcat  $rnaseq_vcf |perl -p -e 's/^/chr/g' | perl -p -e 's/^chr#/#/g' |bgzip > rnaseq.vcf.gz
    zcat somatic_germline.vcf.gz > alltumor.vcf
    tabix rnaseq.vcf.gz
    vcf-merge somatic_germline.vcf.gz rnaseq.vcf.gz |bgzip > allvariants.vcf.gz
    perl $baseDir/vcf2bed.pl alltumor.vcf |cut -f 1,2,3 > tumor.bed
    cat tumor.bed |perl -p -e 's/chr//g' > tumor.nochr.bed
    zgrep '#CHROM' rnaseq.vcf.gz |rev |cut -f 1 |rev  > rnaseqid.txt
    if [[ -z rnaseq.bamreadct.txt ]]
    then
	/project/shared/bicf_workflow_ref/seqprg/bam-readcount/bin/bam-readcount -w 0 -q 0 -b 25 -l tumor.nochr.bed -f ${index_path}/hisat_genome.fa $rbam > rnaseq.bamreadct.txt
	
    fi
    
else
    ln -s somatic_germline.vcf.gz allvariants.vcf.gz
fi
tabix allvariants.vcf.gz
perl $baseDir/integrate_vcfs.pl ${subject} $tumor_id $normal_id $index_path
vcf-sort ${tumor_id}.all.vcf |bedtools intersect -header -a stdin -b ${index_path}/UTSWV2.bed  |uniq |bgzip > ${subject}.utsw.vcf.gz
tabix ${subject}.utsw.vcf.gz
bcftools view -Oz -o ${subject}.vcf.gz -s ${tumor_id}  ${subject}.utsw.vcf.gz
#Makes TumorMutationBurenFile
zgrep -c -v "SS=2" ${subject}.utsw.vcf.gz |awk '{print "Class,TMB\n,"sprintf("%.2f",$1/4.6)}' > ${subject}.tmb.txt

#Convert to HG37
module load crossmap/0.2.5
CrossMap.py vcf ${index_path}/../hg38ToHg19.over.chain.gz ${subject}.utsw.vcf.gz /project/apps_database/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa ${subject}.PASS.hg19.vcf

perl $baseDir/philips_excel.pl ${subject}.PASS.hg19.vcf $rnaseq_fpkm

python $baseDir/../IntellispaceDemographics/gatherdemographics.py -i $subject -u phg_workflow -p $password -o ${subject}.xml
perl /project/PHG/PHG_Clinical/validation/analysis/compareTumorNormal.pl $prefix\.utsw.vcf.gz
