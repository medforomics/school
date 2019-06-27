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
while getopts :r:n:a:p:t:m:n:v:s:i:x:c:d:e:b:f:h opt
do
    case $opt in
        r) index_path=$OPTARG;;
        n) caseID=$OPTARG;;
	a) archive=$OPTARG;;
        p) subject=$OPTARG;;
        t) tumor_id=$OPTARG;;
        n) normal_id=$OPTARG;;
        v) tumor_vcf=$OPTARG;;
        s) somatic_vcf=$OPTARG;;
	m) merged_vcf=$OPTARG;;
	i) itd_vcf=$OPTARG;;
	d) cnv_answer=$OPTARG;;
        x) rnaseq_vcf=$OPTARG;;
        c) rnaseq_ntct=$OPTARG;;
	f) rnaseq_fpkm=$OPTARG;;
	e) rnaseq_bam=$OPTARG;;
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
targetbed="${index_path}/clinseq_prj/UTSWV2_2.panelplus.bed"
fi

baseDir="`dirname \"$0\"`"
vepdir='/project/shared/bicf_workflow_ref/vcf2maf'

module load bedtools/2.26.0 samtools/gcc/1.8 htslib/gcc/1.8 vcftools/0.1.14 snpeff/4.3q

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
    if [[ -f "somatic_${dna_runid}/${subject}_${dna_runid}.somatic.vcf.gz" ]]
    then
	somaticvcf="somatic_${dna_runid}/${subject}_${dna_runid}.somatic.vcf.gz"
    fi
    if [[ -f "${subject}_${dna_runid}.germline.vcf.gz" ]] 
    then
	tumor_vcf="${subject}_${dna_runid}.germline.vcf.gz"
    fi
    if [[ -d  "${subject}_${dna_runid}.dna.vcf.gz" ]]
    then
	merged_vcf="${subject}_${dna_runid}.dna.vcf.gz"
    fi
    if [[ -f "${subject}.pindel_tandemdup.pass.vcf.gz" ]]
    then
	itd_vcf="${subject}.pindel_tandemdup.pass.vcf.gz"
    elif [[ -f "dna_${dna_runid}/${subject}.pindel_tandemdup.pass.vcf.gz" ]]
    then
	itd_vcf="dna_${dna_runid}/${subject}.pindel_tandemdup.pass.vcf.gz"
    fi
    cnv_answer="$tumor_id/$tumor_id.cnv.answer.txt"
    
    if [[ -f "${subject}_${rna_runid}*rna.vcf.gz" ]]
    then
	if [[ -f "${subject}_${rna_runid}.germline.rna.vcf.gz" ]]
	then
	    rnaseq_vcf="${subject}_${rna_runid}.germline.rna.vcf.gz"
	elif [[ -f "${subject}_${rna_runid}.rna.vcf.gz" ]]
	then
	    rnaseq_vcf="${subject}_${rna_runid}.rna.vcf.gz"
	fi
	rnaseq_ntct="${rnaseq_id}/${rnaseq_id}.bamreadct.txt"
	rnaseq_fpkm="${rnaseq_id}/${rnaseq_id}.fpkm.txt"
	rnaseq_bam="${rnaseq_id}/${rnaseq_id}.bam"
    fi
fi

    #Merge VCFs if there are 2
if [[ -f $somatic_vcf && -f $tumor_vcf ]]
then
    tabix -f $somatic_vcf
    tabix -f $tumor_vcf
    bcftools annotate -Ov -a $tumor_vcf -o somatic.only.vcf --columns CHROM,POS,CallSet $somatic_vcf
    vcf-shuffle-cols -t somatic.only.vcf $tumor_vcf |bgzip > tumor.vcf.gz
    bgzip -f somatic.only.vcf
    vcf-concat somatic.only.vcf.gz tumor.vcf.gz |vcf-sort |uniq | bgzip > somatic_germline.vcf.gz
elif [[ -f $merged_vcf ]]
then
    cp $merged_vcf somatic_germline.vcf.gz
elif [[ -f $tumor_vcf ]]
then   
    cp $tumor_vcf somatic_germline.vcf.gz
else
    cp $somatic_vcf somatic_germline.vcf.gz
fi
tabix -f somatic_germline.vcf.gz

#Merge ITD with CNV File

if [[ -f $itd_vcf ]]
then
    perl $baseDir/itdvcf2cnv.pl $tumor_id $itd_vcf
    cat dupcnv.txt >> $cnv_answer
fi

#Filter VCF File
icommand="perl $baseDir/integrate_vcfs.pl -s ${subject} -t $tumor_id -r $index_path"
if [[ -n $normal_id ]]
then 
    icommand+=" -n $normal_id"
fi
if [[ -f $rnaseq_vcf ]]
then
    icommand+=" -v $rnaseq_vcf -c $rnaseq_ntct"
fi
$icommand
vcf-sort ${subject}.all.vcf | bedtools intersect -header -a stdin -b $targetbed | uniq | bgzip > ${subject}.vcf.gz
bgzip -f ${subject}.pass.vcf
tabix -f ${subject}.vcf.gz
tabix -f ${subject}.pass.vcf.gz

#Identify functional splice sites
if [[ -f ${rnaseq_bam} ]]
then
    samtools index -@ 4 ${rnaseq_bam}
    /project/shared/bicf_workflow_ref/seqprg/bin/regtools cis-splice-effects identify ${subject}.vcf.gz ${rnaseq_bam} ${index_path}/genome.fa ${index_path}/gencode.gtf -o ${subject}.splicevariants.txt -v ${subject}.splicevariants.vcf -s 0 -e 5 -i 5
fi

#Makes TumorMutationBurenFile

bedtools intersect -header -a ${subject}.pass.vcf.gz -b $targetbed |uniq |bgzip > ${subject}.utswpass.vcf.gz

#Calculate TMB
targetsize=`awk '{sum+=$3-$2} END {print sum/1000000}' $targetbed`
if [[ -n $normal_id ]]
then
zgrep "#\|SS=2" ${subject}.utswpass.vcf.gz > ${subject}.utswpass.somatic.vcf
zgrep -c -v "#" ${subject}.utswpass.somatic.vcf.gz | awk -v tsize="$targetsize" '{print "Class,TMB\n,"sprintf("%.2f",$1/tsize)}' > ${subject}.TMB.csv
perl $baseDir/compareTumorNormal.pl ${subject}.utswpass.vcf.gz > ${subject}.concordance.txt
else
    echo -e "Class,TMB\n,0.00" > ${subject}.TMB.csv
fi

if [[ -a $archive ]]
then
    bash $baseDir/syncCase2Azure.sh ${subject}
fi
