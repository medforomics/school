#!/bin/bash
#union.sh

usage() {
  echo "-h Help documentation for gatkrunner.sh"
  echo "-r  --Reference Genome: GRCh38 or GRCm38"
  echo "-p  --Prefix for output file name"
  echo "-n  --NuCLIA CaseID"
  echo "-t  --TumorID"
  echo "-v  --tumor vcf"
  echo "-s  --somatic vcf"
  echo "-i  --indel vcf"
  echo "-x  --rnaseq vcf"
  echo "-c  --rnaseq read ct"
  echo "-f  --rnaseq fpkm"
  echo "-b  --targetbed"
  echo "-a  --archive"
  echo "Example: bash unify_case.sh -p prefix -r /path/GRCh38"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :r:n:a:p:t:m:g:v:s:k:y:j:i:x:c:d:e:b:f:h opt
do
    case $opt in
        r) index_path=$OPTARG;;
        y) rna_ref_path=$OPTARG;;
        n) caseID=$OPTARG;;
	a) archive=$OPTARG;;
        p) subject=$OPTARG;;
        t) tumor_id=$OPTARG;;
        g) normal_id=$OPTARG;;
        v) tumor_vcf=$OPTARG;;
        s) somatic_vcf=$OPTARG;;
	m) merged_vcf=$OPTARG;;
	i) itd_vcf=$OPTARG;;
	d) cnv_answer=$OPTARG;;
        x) rnaseq_vcf=$OPTARG;;
        c) rnaseq_ntct=$OPTARG;;
	f) rnaseq_fpkm=$OPTARG;;
	e) rnaseq_bam=$OPTARG;;
	j) rnaseq_translocation=$OPTARG;;
 	b) targetbed=$OPTARG;;
	k) nucliatoken=$OPTARG;;
	h) usage;;
    esac
done

shift $(($OPTIND -1))

module load bedtools/2.26.0 samtools/gcc/1.8 bcftools/gcc/1.8 htslib/gcc/1.8 vcftools/0.1.14 snpeff/4.3q
if [[ -z $targetbed ]]
then
targetbed="${index_path}/../panels/UTSW_V4_heme/targetpanel.bed"
fi

baseDir="`dirname \"$0\"`"
vepdir='/project/shared/bicf_workflow_ref/vcf2maf'

dna_ref_path=$index_path
if [[ -z $rna_ref_path ]]
then
    rna_ref_path="${index_path}/../rnaref"
fi

gfopt=''
svcalls=''
if [[ -n ${caseID} ]]
then
    if [[ ! -f ${caseID}.json ]]
    then
	curl -o ${caseID}.json "https://nuclia.biohpc.swmed.edu/getSampleAndRunIds?token=$nucliatoken&projectId=$caseID"
    fi
    jsonout=$(perl ${baseDir}/parse_getsampleJson.pl ${caseID}.json 2>&1)
    myarray=($jsonout)
    subject=${myarray[1]}
    tumor_id=${myarray[3]}
    dna_runid=${myarray[2]}
    rna_runid=${myarray[6]}
    normal_id=${myarray[5]}
    rnaseq_id=${myarray[7]}

    if [[ -f  "${subject}_${dna_runid}.dna.vcf.gz" ]]
    then
	merged_vcf="${subject}_${dna_runid}.dna.vcf.gz"
    fi
    if [[ -n `ls ${subject}_${rna_runid}*rna.vcf.gz` ]]
    then
	rnaseq_vcf="${subject}_${rna_runid}.rna.vcf.gz"
	rnaseq_ntct="${rnaseq_id}/${rnaseq_id}.bamreadct.txt"
	rnaseq_fpkm="${rnaseq_id}/${rnaseq_id}.fpkm.txt"
	rnaseq_bam="${rnaseq_id}/${rnaseq_id}.bam"
	rnaseq_translocation="${rnaseq_id}/${rnaseq_id}.translocations.answer.txt"
	gfopt="${gfopt} -f ${rnaseq_translocation}"
    fi
    cnv_answer="$tumor_id/$tumor_id.cnv.answer.txt"
fi

echo $subject $tumor_id $tumor_vcf $somatic_vcf $merged_vcf 
if [[ -z $subject ]] || [[ -z $index_path ]]; then
    usage
fi 

#Merge VCFs if there are 2

cp $merged_vcf somatic_germline.vcf.gz
tabix -f somatic_germline.vcf.gz

#Identify functional splice sites

if [[ -f ${rnaseq_bam} ]]
then
    samtools index -@ 4 ${rnaseq_bam}
    /project/shared/bicf_workflow_ref/seqprg/bin/regtools cis-splice-effects identify somatic_germline.vcf.gz ${rnaseq_bam} ${index_path}/genome.fa ${index_path}/gencode.gtf -o ${subject}.splicevariants.txt -v ${subject}.splicevariants.vcf -s 0 -e 5 -i 5
fi
#Filter VCF File
icommand="perl $baseDir/integrate_vcfs.pl -s ${subject} -t $tumor_id -r $index_path"
if [[ -n $normal_id ]]
then 
    icommand+=" -n $normal_id"
fi
if [[ -f $rnaseq_vcf ]]
then
    icommand+=" -v $rnaseq_vcf -c $rnaseq_ntct -g ${subject}.splicevariants.vcf"
fi
$icommand

vcf-sort ${subject}.all.vcf | bedtools intersect -header -a stdin -b $targetbed | uniq | bgzip > ${subject}.set1.vcf.gz
bgzip -f ${subject}.pass.vcf
tabix -f ${subject}.pass.vcf.gz

if [[ -f "${subject}.pindel_tandemdup.vcf.gz" ]]
then
    itdpindel_vcf="${subject}.pindel_tandemdup.vcf.gz"
    gfopt="${gfopt} -i ${subject}.pindel.genefusion.txt"
    if [[ -n $rnaseq_id ]]
    then
	cp ${subject}.pindel.vcf.gz ${subject}.pindel.ori.vcf.gz
	perl ${baseDir}/add_blank_sample_vcf.pl -r $rnaseq_id -v ${subject}.pindel.ori.vcf.gz -o ${subject}.pindel.vcf
	bgzip -f ${subject}.pindel.vcf
    fi
    vcf-shuffle-cols -t ${subject}.set1.vcf.gz ${subject}.pindel.vcf.gz | bgzip > pindel.vcf.gz
    svcalls="${svcalls} pindel.vcf.gz"
elif [[ -f "dna_${dna_runid}/${subject}.pindel_tandemdup.vcf.gz" ]]
then
    itdpindel_vcf="dna_${dna_runid}/${subject}.pindel_tandemdup.vcf.gz"
    gfopt="${gfopt} -i dna_${dna_runid}/${subject}.pindel.genefusion.txt"
    if [[ -n $rnaseq_id ]]
    then
	cp dna_${dna_runid}/${subject}.pindel.vcf.gz dna_${dna_runid}/${subject}.pindel.ori.vcf.gz
	perl ${baseDir}/add_blank_sample_vcf.pl -r $rnaseq_id -v dna_${dna_runid}/${subject}.pindel.ori.vcf.gz -o dna_${dna_runid}/${subject}.pindel.vcf
	bgzip -f dna_${dna_runid}/${subject}.pindel.vcf
    fi
    vcf-shuffle-cols -t ${subject}.set1.vcf.gz dna_${dna_runid}/${subject}.pindel.vcf.gz | bgzip > pindel.vcf.gz
    svcalls="${svcalls} pindel.vcf.gz"
fi
if [[ -f "${tumor_id}/${tumor_id}.itdseek_tandemdup.vcf.gz" ]]
then
    itdseeker_vcf="${tumor_id}/${tumor_id}.itdseek_tandemdup.vcf.gz"
fi
if [[ -f "dna_${dna_runid}/${subject}.delly.vcf.gz" ]]
then
    gfopt="${gfopt} -d dna_${dna_runid}/${subject}.delly.genefusion.txt"
    if [[ -n $rnaseq_id ]]
    then
	cp dna_${dna_runid}/${subject}.delly.vcf.gz dna_${dna_runid}/${subject}.delly.ori.vcf.gz
	perl ${baseDir}/add_blank_sample_vcf.pl -r $rnaseq_id -v dna_${dna_runid}/${subject}.delly.ori.vcf.gz -o dna_${dna_runid}/${subject}.delly.vcf
	bgzip -f dna_${dna_runid}/${subject}.delly.vcf
    fi
    vcf-shuffle-cols -t ${subject}.set1.vcf.gz dna_${dna_runid}/${subject}.delly.vcf.gz | bgzip > delly.vcf.gz
    svcalls="${svcalls} delly.vcf.gz"
fi
if [[ -f "dna_${dna_runid}/${subject}.svaba.vcf.gz" ]]
then
    gfopt="${gfopt} -s dna_${dna_runid}/${subject}.svaba.genefusion.txt"
    if [[ -n $rnaseq_id ]]
    then
	cp dna_${dna_runid}/${subject}.svaba.vcf.gz dna_${dna_runid}/${subject}.svaba.ori.vcf.gz
	perl ${baseDir}/add_blank_sample_vcf.pl -r $rnaseq_id -v dna_${dna_runid}/${subject}.svaba.ori.vcf.gz -o dna_${dna_runid}/${subject}.svaba.vcf
	bgzip -f dna_${dna_runid}/${subject}.svaba.vcf
    fi
    vcf-shuffle-cols -t ${subject}.set1.vcf.gz dna_${dna_runid}/${subject}.svaba.vcf.gz | bgzip > svaba.vcf.gz
    svcalls="${svcalls} svaba.vcf.gz"
fi

vcf-concat  ${subject}.set1.vcf.gz $svcalls |vcf-sort |bgzip > ${subject}.vcf.gz
tabix -f ${subject}.vcf.gz

#merge gene fusion files
perl $baseDir/merge_genefusions.pl -p ${subject} -t ${tumor_id} -r ${rna_ref_path} $gfopt


#Merge ITD with CNV File
if [[ -f $itdpindel_vcf ]] && [[ -f $itdseeker_vcf ]]
then
    tabix -f $itdpindel_vcf
    tabix -f $itdseeker_vcf
    bcftools view -Oz -o itdpindel.vcf.gz -s $tumor_id $itdpindel_vcf
    tabix itdpindel.vcf.gz
    bcftools concat -Oz -o itd.vcf.gz itdpindel.vcf.gz  $itdseeker_vcf
    perl $baseDir/itdvcf2cnv.pl $tumor_id itd.vcf.gz
    cat ${tumor_id}.dupcnv.txt >> $cnv_answer
elif [[ -f $itdpindel_vcf ]]
then
    perl $baseDir/itdvcf2cnv.pl $tumor_id $itdpindel_vcf
    cat ${tumor_id}.dupcnv.txt >> $cnv_answer
#fi
elif [[ -f $itdseeker_vcf ]]
then
    perl $baseDir/itdvcf2cnv.pl $tumor_id $itdseeker_vcf
    cat ${tumor_id}.dupcnv.txt >> $cnv_answer
fi

#Makes TumorMutationBurenFile

bedtools intersect -header -a ${subject}.pass.vcf.gz -b $targetbed |uniq |bgzip > ${subject}.utswpass.vcf.gz
targetsize=`awk '{sum+=$3-$2} END {print sum/1000000}' $targetbed`

if [[ -n $normal_id ]] && [[ $normal_id != 'NA' ]]
then
zgrep "#\|SS=2" ${subject}.utswpass.vcf.gz > ${subject}.utswpass.somatic.vcf
grep -c -v "#" ${subject}.utswpass.somatic.vcf | awk -v tsize="$targetsize" '{print "Class,TMB\n,"sprintf("%.2f",$1/tsize)}' > ${subject}.TMB.csv
perl $baseDir/compareTumorNormal.pl ${subject}.utswpass.vcf.gz > ${subject}.concordance.txt
perl ${vepdir}/vcf2maf.pl --input ${subject}.utswpass.somatic.vcf --output ${caseID}.maf --species homo_sapiens --ncbi-build GRCh38 --ref-fasta ${vepdir}/.vep/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa --filter-vcf ${vepdir}/.vep/homo_sapiens/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz --cache-version 91 --vep-path ${vepdir}/variant_effect_predictor --tumor-id $tumor_id --normal-id $normal_id --custom-enst ${vepdir}/data/isoform_overrides_uniprot --custom-enst ${vepdir}/data/isoform_overrides_at_mskcc --maf-center http://www.utsouthwestern.edu/sites/genomics-molecular-pathology/ --vep-data ${vepdir}/.vep
else
    echo -e "Class,TMB\n,0.00" > ${subject}.TMB.csv
fi

#RUN VCF2MAF
perl ${vepdir}/vcf2maf.pl --input ${subject}.all.vcf --output ${caseID}.maf --species homo_sapiens --ncbi-build GRCh38 --ref-fasta ${vepdir}/.vep/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa --filter-vcf ${vepdir}/.vep/homo_sapiens/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz --cache-version 91 --vep-path ${vepdir}/variant_effect_predictor --tumor-id $tumor_id --normal-id $normal_id --custom-enst ${vepdir}/data/isoform_overrides_uniprot --custom-enst ${vepdir}/data/isoform_overrides_at_mskcc --maf-center http://www.utsouthwestern.edu/sites/genomics-molecular-pathology/ --vep-data ${vepdir}/.vep

if [[ -n $archive ]]
then
    bash $baseDir/syncCase2Azure.sh ${subject} $nucliatoken
fi
