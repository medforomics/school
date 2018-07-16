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
  echo "-c  --rnaseq read ct"
  echo "-f  --rnaseq fpkm"
  echo "-b  --targetbed"
  echo "Example: bash unify_case.sh -p prefix -r /path/GRCh38"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :r:p:t:n:v:s:x:c:b:f:h opt
do
    case $opt in
        r) index_path=$OPTARG;;
        p) subject=$OPTARG;;
        t) tumor_id=$OPTARG;;
        n) normal_id=$OPTARG;;
        v) tumor_vcf=$OPTARG;;
        s) somatic_vcf=$OPTARG;;
        x) rnaseq_vcf=$OPTARG;;
        c) rnaseq_ntct=$OPTARG;;
	f) rnaseq_fpkm=$OPTARG;;
 	b) targetbed=$OPTARG;;
	h) usage;;
    esac
done

shift $(($OPTIND -1))

if [[ -z $tumor_vcf ]] || [[ -z $subject ]] || [[ -z $index_path ]]; then
    usage
fi 
if [[ -z $targetbed ]]
then
targetbed="${index_path}/clinseq_prj/UTSWV2.bed"
fi

baseDir="`dirname \"$0\"`"
vepdir='/project/shared/bicf_workflow_ref/vcf2maf'

module load bedtools/2.26.0 samtools/1.6 vcftools/0.1.14
if [[ -a $somatic_vcf ]] 
then
    tabix -f $somatic_vcf
    tabix -f $tumor_vcf
    bcftools annotate -Ov -a $tumor_vcf -o somatic.only.vcf --columns CHROM,POS,CallSet $somatic_vcf
    vcf-shuffle-cols -t somatic.only.vcf $tumor_vcf |bgzip > tumor.vcf.gz
    bgzip -f somatic.only.vcf
    vcf-concat somatic.only.vcf.gz tumor.vcf.gz |vcf-sort |uniq | bgzip > somatic_germline.vcf.gz
else
    ln -s $tumor_vcf somatic_germline.vcf.gz
fi

tabix -f somatic_germline.vcf.gz

perl $baseDir/integrate_vcfs.pl -s ${subject} -t $tumor_id -n $normal_id -r $index_path -v $rnaseq_vcf -c $rnaseq_ntct
vcf-sort ${subject}.all.vcf | bedtools intersect -header -a stdin -b $targetbed | uniq | bgzip > ${subject}.vcf.gz
bgzip -f ${subject}.pass.vcf
tabix -f ${subject}.vcf.gz
tabix -f ${subject}.pass.vcf.gz

#Makes TumorMutationBurenFile

bedtools intersect -header -a ${subject}.pass.vcf.gz -b $targetbed |uniq |bgzip > ${subject}.utswpass.vcf.gz

targetsize=`awk '{sum+=$3-$2} END {print sum/1000000}' $targetbed`
zgrep -c "SS=2" ${subject}.utswpass.vcf.gz | awk -v tsize="$targetsize" '{print "Class,TMB\n,"sprintf("%.2f",$1/tsize)}' > ${subject}.TMB.csv
perl $baseDir/compareTumorNormal.pl ${subject}.utswpass.vcf.gz > ${subject}.concordance.txt

#Convert to HG37
module load crossmap/0.2.5
CrossMap.py vcf ${index_path}/../hg38ToHg19.over.chain.gz ${subject}.utswpass.vcf.gz /project/apps_database/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa ${subject}.PASS.hg19.vcf
cat ${subject}.PASS.hg19.vcf |perl -p -e 's/^chr//g' > ${subject}.formaf.vcf
perl $baseDir/philips_excel.pl ${subject}.PASS.hg19.vcf $rnaseq_fpkm

zcat ${subject}.pass.vcf.gz |perl -p -e 's/^chr//g' > ${subject}.formaf.vcf
perl ${vepdir}/vcf2maf.pl --input ${subject}.formaf.vcf --output ${subject}.maf --species homo_sapiens --ncbi-build GRCh37 --ref-fasta ${vepdir}/.vep/homo_sapiens/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa --filter-vcf ${vepdir}/.vep/homo_sapiens/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz --cache-version 91 --vep-path ${vepdir}/variant_effect_predictor --tumor-id $tumor_id --normal-id $normal_id --custom-enst ${vepdir}/data/isoform_overrides_uniprot.txt --custom-enst ${vepdir}/data/isoform_overrides_at_mskcc.txt --maf-center http://www.utsouthwestern.edu/sites/genomics-molecular-pathology/ --vep-data ${vepdir}/.vep

rsync -rlptgoD ${subject}.vcf.gz ${subject}.TMB.csv /project/PHG/PHG_BarTender/bioinformatics/seqanalysis/${subject}

python $baseDir/../IntellispaceDemographics/gatherdemographics.py -i $subject -u phg_workflow -p $password -o ${subject}.xml

#sequencestats
#exoncoverage
