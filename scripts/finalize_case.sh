#!/bin/bash
#union.sh

usage() {
  echo "-h Help documentation for gatkrunner.sh"
  echo "-r  --Reference Genome: GRCh38 or GRCm38"
  echo "-p  --Prefix for output file name"
  echo "-n  --NuCLIA CaseID"
  echo "-b  --targetbed"
  echo "-a  --archive"
  echo "Example: bash unify_case.sh -p prefix -r /path/GRCh38"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :n:r:y:a:p:d:t:g:b:k:l:h opt
do
    case $opt in
        n) caseID=$OPTARG;;
        r) index_path=$OPTARG;;
	p) seqrunid=$OPTARG;;
	d) casedir=$OPTARG;;
        y) rna_ref_path=$OPTARG;;
	a) archive=$OPTARG;;
        t) tumor_id=$OPTARG;;
        g) normal_id=$OPTARG;;
 	b) targetbed=$OPTARG;;
	k) nucliatoken=$OPTARG;;
	l) testdir=$OPTARG;;
	h) usage;;
    esac
done

shift $(($OPTIND -1))

module load bedtools/2.26.0 samtools/gcc/1.8 bcftools/gcc/1.8 htslib/gcc/1.8 vcftools/0.1.14 snpeff/4.3q R/3.6.1-gccmkl

if [[ -z $index_path ]]
then
    index_path=/project/shared/bicf_workflow_ref/human/grch38_cloud/dnaref
fi
dna_ref_path=$index_path
if [[ -z $rna_ref_path ]]
then
    rna_ref_path="${index_path}/../rnaref"
fi

baseDir="`dirname \"$0\"`"
cd $casedir
echo "*****Set File Variables******"
normal_id='NA'
gfout=''
for bamf in *naout/*.bam
do
    fname=$(basename $bamf)
    if [[ $bamf =~ 'T_DNA' ]] && [[ $bamf =~ 'consensus' ]]
    then
	tumor_id="${fname%.consensus.bam}"
	tumorbam=$bamf
    elif [[ $bamf =~ 'T_RNA' ]]
    then
	rnaseq_id="${fname%.bam}"
	rnavcf="rnavcf/${caseID}.fb.vcf.gz"
	rnaseq_ntct="rnaout/${rnaseq_id}.bamreadcount.txt.gz"
	rnaseq_fpkm="rnaout/${rnaseq_id}.fpkm.txt"
	rnaseq_bam="rnaout/${rnaseq_id}.bam"
	rnaseq_translocation="rnaout/${rnaseq_id}.translocations.answer.txt"
	gfopt="${gfopt} -f ${rnaseq_translocation}"
    elif [[ $bamf =~ 'N_DNA' ]]
    then
	normal_id=$sid
    fi
done

dnavcf="dnavcf/${caseID}.union.vcf.gz"
cnv_answer="dnaout/$tumor_id.cnv.answer.txt"

echo $caseID $tumor_id $dnavcf 
echo "*****Load Data into NuCLIA******"
if [[ -z $caseID ]] || [[ -z $index_path ]] || [[ -z $dnavcf ]]
then
    usage
fi 

if [[ -z $testdir ]]
then
    ssh answerbe@198.215.54.71 "mkdir -p /swnas/qc_nuclia/seqanalysis/${seqrunid}"
    if [[ $caseID =~ 'GM12878' ]]
    then
	bash ${baseDir}/scripts/snsp.sh -p $tumor_id -t $targetbed -b $tumorbam -v $dnavcf
    fi
    perl $baseDir/create_propertiesfile.pl $casedir
    rsync -rlptgoD --exclude="*fastq.gz*" --exclude "*work*" --exclude="*bam*" $casedir answerbe@198.215.54.71:/swnas/qc_nuclia/seqanalysis/${seqrunid}
    for j in *properties ; do
	pfile=$(basename $j)
	curl "http://nuclia.biohpc.swmed.edu:8080/NuCLIAVault/addPipelineResultsWithProp?token=${nucliatoken}&propFilePath=/swnas/qc_nuclia/seqanalysis/${seqrunid}/${pfile}"
    done
fi

echo "*****Concat Viral RNA Results******"
if [[ -f "dnaout/${tumor_id}.viral.seqstats.txt" ]]
then
    perl $baseDir/parse_viral_idxstats.pl dnaout/*.viral.seqstats.txt
    mv viral_results.txt ${caseID}.viral_results.txt
fi

vepdir='/project/shared/bicf_workflow_ref/vcf2maf'
svcalls=''

echo "*****Creating a VCF File******"
#Identify functional splice sites

icommand="perl $baseDir/integrate_vcfs.pl -f $dnavcf -s ${caseID} -t $tumor_id -r $index_path"
if [[ -n $normal_id ]]
then 
    icommand+=" -n $normal_id"
fi
if [[ -f $rnavcf ]]
then
    samtools index -@ 4 ${rnaseq_bam}
    /project/shared/bicf_workflow_ref/seqprg/bin/regtools cis-splice-effects identify ${dnavcf} ${rnaseq_bam} ${index_path}/genome.fa ${index_path}/gencode.gtf -o ${caseID}.splicevariants.txt -v ${caseID}.splicevariants.vcf -s 0 -e 5 -i 5
    icommand+=" -v $rnavcf -c $rnaseq_ntct -g ${caseID}.splicevariants.vcf"
else
    unset rnaseq_id
fi
$icommand

vcf-sort ${caseID}.all.vcf | bedtools intersect -header -a stdin -b $targetbed | uniq | bgzip > ${caseID}.set1.vcf.gz
bgzip -f ${caseID}.pass.vcf
tabix -f ${caseID}.pass.vcf.gz

echo "*****Evaluating ITDs******"
if [[ -f "dnacallset/${caseID}.pindel_tandemdup.vcf.gz" ]]
then
    itdpindel_vcf="dnacallset/${caseID}.pindel_tandemdup.vcf.gz"
    gfopt="${gfopt} -i dnacallset/${caseID}.pindel.genefusion.txt"
fi
if [[ -f "dnacallset/${tumor_id}.itdseek_tandemdup.vcf.gz" ]]
then
    itdseeker_vcf="dnacallset/${tumor_id}.itdseek_tandemdup.vcf.gz"
fi
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
elif [[ -f $itdseeker_vcf ]]
then
    perl $baseDir/itdvcf2cnv.pl $tumor_id $itdseeker_vcf
    cat ${tumor_id}.dupcnv.txt >> $cnv_answer
fi

echo "*****Adding SVs******"
svcallers="delly svaba.sv"
for c in $svcallers
do	 
    gfopt="${gfopt} dnacallset/${caseID}.${c}.genefusion.txt"
    if [[ -n $rnaseq_id ]]
    then
	perl ${baseDir}/add_blank_sample_vcf.pl -t $tumor_id -n $normal_id -r $rnaseq_id -v dnacallset/${caseID}.${c}.vcf.gz -o dnacallset/${caseID}.${c}.final.vcf
    elif [[ $normal_id == 'NA' ]]
    then
	perl ${baseDir}/add_blank_sample_vcf.pl -t $tumor_id -n $normal_id -v dnacallset/${caseID}.${c}.vcf.gz -o dnacallset/${caseID}.${c}.final.vcf
	bgzip -f dnacallset/${caseID}.${c}.final.vcf
    else
    	zcat dnacallset/${caseID}.${c}.vcf.gz | perl -pe 's/${caseID}.tumor/${tumor_id}/g' | perl -pe 's/${caseID}.normal/${normal_id}/g' > dnacallset/${caseID}.${c}.final.vcf.gz
    fi
    bgzip -f dnacallset/${caseID}.${c}.final.vcf
    tabix dnacallset/${caseID}.${c}.final.vcf.gz
    vcf-shuffle-cols -t ${caseID}.set1.vcf.gz dnacallset/${caseID}.${c}.final.vcf.gz | bgzip > ${c}.vcf.gz
    svcalls="${svcalls} ${c}.vcf.gz"
done

#Merge SV to SNVs/Indel
vcf-concat  ${caseID}.set1.vcf.gz $svcalls |vcf-sort |bgzip > ${caseID}.vcf.gz
tabix -f ${caseID}.vcf.gz

#merge gene fusion files
perl $baseDir/merge_genefusions.pl -p ${caseID} -t ${tumor_id} -r ${rna_ref_path} $gfopt

echo "*****Calculate TMB and MSI******"
#Makes TumorMutationBurenFile

bedtools intersect -header -a ${caseID}.pass.vcf.gz -b $targetbed |uniq |bgzip > ${caseID}.utswpass.vcf.gz
targetsize=`awk '{sum+=$3-$2} END {print sum/1000000}' $targetbed`

if [[ -n $normal_id ]] && [[ $normal_id != 'NA' ]]
then
    zgrep "#\|SS=2" ${caseID}.utswpass.vcf.gz > ${caseID}.utswpass.somatic.vcf
    tmbval=`grep -c -v "#" ${caseID}.utswpass.somatic.vcf | awk -v tsize="$targetsize" '{print sprintf("%.2f",$1/tsize)}'`
    tmbclass='Low-TMB'
    if [[ $(echo "$tmbval >= 5" |bc -l) == 1 ]]
    then
	tmbclass='Medium-TMB' 
	java -jar $SNPEFF_HOME/SnpSift.jar filter "(TYPE=='snp')" ${caseID}.utswpass.somatic.vcf > ${caseID}.snps.vcf
	Rscript ${baseDir}/run_mutsig.R
	mv mutational_signature.txt ${caseID}.mutational_signature.txt
	mv mutational_signature.png ${caseID}.mutational_signature.png
    fi
    if [[ $(echo "$tmbval >= 10" |bc -l) == 1 ]]
    then
	tmbclass='High-TMB' 
    fi
    perl $baseDir/compareTumorNormal.pl ${caseID}.utswpass.vcf.gz > ${caseID}.concordance.txt
    vep_opt="--normal-id $normal_id"
else
    tmbval=0.00
    tmbclass=UNK
fi

msiclass='MSS'
msival=`grep -v Total dnacallset/${caseID}.msi |cut -f 3`
if [[ $(echo "$msival >= 10" |bc -l) == 1 ]]
then
    msiclass='MSI'
fi

echo -e "Metric,Value,Class\nTMB,${tmbval},${tmbclass}\nMSI,${msival},${msiclass}" > ${caseID}.TMB.csv


#RUN VCF2MAF
perl ${vepdir}/vcf2maf.pl --input ${caseID}.all.vcf --output ${caseID}.maf --species homo_sapiens --ncbi-build GRCh38 --ref-fasta ${vepdir}/.vep/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa --filter-vcf ${vepdir}/.vep/homo_sapiens/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz --cache-version 91 --vep-path ${vepdir}/variant_effect_predictor --tumor-id $tumor_id $vep_opt --custom-enst ${vepdir}/data/isoform_overrides_uniprot --custom-enst ${vepdir}/data/isoform_overrides_at_mskcc --maf-center http://www.utsouthwestern.edu/sites/genomics-molecular-pathology/ --vep-data ${vepdir}/.vep

if [[ -n $archive ]]
then
    if [[ -n $testdir ]]
    then
	wkdir=$testdir
	testing=1
    else
	wkdir=/archive/PHG/PHG_Clinical
    fi
    mkdir -p ${wkdir}/cases/${caseID}
    if [[ -f ${caseID}.vcf.gz ]]
    then
	cp ${caseID}.vcf.gz ${wkdir}/cases/${caseID}
	cp ${caseID}.maf ${wkdir}/cases/${caseID}
	cp ${caseID}.viral_results.txt ${wkdir}/cases/${caseID}
	cp ${caseID}.TMB.csv ${wkdir}/cases/${caseID}
    fi
    if [[ -f ${caseID}.mutational_signature.txt ]]
    then
	cp ${caseID}.mutational_signature.txt ${wkdir}/cases/${caseID}
	cp ${caseID}.mutational_signature.png ${wkdir}/cases/${caseID}
    fi
    if [[ -f  dnacallset/${caseID}.cnv.answer.txt ]]
    then
	mkdir -p ${wkdir}/cases/${caseID}/${tumor_id}
	cp dnacallset/${caseID}.answerplot.cnr ${wkdir}/cases/${caseID}/${tumor_id}/${tumor_id}.answerplot.cnr
	cp dnacallset/${caseID}.answerplot.cns ${wkdir}/cases/${caseID}/${tumor_id}/${tumor_id}.answerplot.cns
	cp dnacallset/${caseID}.cnv.answer.txt ${wkdir}/cases/${caseID}/${tumor_id}/${tumor_id}.cnv.answer.txt
	cp dnacallset/${caseID}.ballelefreq.txt ${wkdir}/cases/${caseID}/${tumor_id}/${tumor_id}.ballelefreq.txt
    fi
    if [[ -f ${caseID}.translocations.answer.txt ]]
    then
	cp ${caseID}.translocations.answer.txt ${wkdir}/cases/${caseID}
    fi
    if [[ -f dnaout/${tumor_id}.consensus.bam ]]
    then
	cp dnaout/${tumor_id}.consensus.bam ${wkdir}/cases/${caseID}/${tumor_id}/${tumor_id}.bam
	samtools index -@ 4 ${wkdir}/cases/${caseID}/${tumor_id}/${tumor_id}.bam
    fi
    if [[ -f dnaout/${normal_id}.consensus.bam ]]
    then
	mkdir -p ${wkdir}/cases/${caseID}/${normal_id}
	cp ${caseID}.concordance.txt ${wkdir}/cases/${caseID}
	cp ${caseID}.utswpass.somatic.vcf* ${wkdir}/cases/${caseID}
	if [[ -f ${wkdir}/cases/${caseID}/${caseID}.utswpass.somatic.vcf.gz ]]
	then
	    gunzip ${wkdir}/cases/${caseID}/${caseID}.utswpass.somatic.vcf.gz
	fi
	cp dnaout/${normal_id}.consensus.bam ${wkdir}/cases/${caseID}/${normal_id}/${normal_id}.bam
	samtools index -@ 4 ${wkdir}/cases/${caseID}/${normal_id}/${normal_id}.bam
    fi	   
    if [[ -f rnaout/${rnaseq_id}.bam ]]
    then
	mkdir -p ${wkdir}/cases/${caseID}/${rnaseq_id}
	cp rnaout/${rnaseq_id}.bam ${wkdir}/cases/${caseID}/${rnaseq_id}
	samtools index -@ 4 ${wkdir}/cases/${caseID}/${rnaseq_id}/${rnaseq_id}.bam
	cp rnaout/${rnaseq_id}.cts ${wkdir}/cases/${caseID}/${rnaseq_id}
	cp rnaout/${rnaseq_id}.fpkm.txt ${wkdir}/cases/${caseID}/${rnaseq_id}
	if [[ ! -f ${caseID}.translocations.answer.txt ]]
	then
            cp rnaout/${rnaseq_id}.translocations.answer.txt ${wkdir}/cases/${caseID}/${rnaseq_id}
	fi
    fi
    if [[ -z $testing ]]
    then
	mv $casedir $wkdir/toarchive/caseDirs/${caseID}
	rsync -avz $wkdir/cases/${caseID} answerbe@198.215.54.71:/swnas/cases
	source ${baseDir}/../azure_credentials
	az storage blob upload-batch -d cases -s $wkdir/cases/${caseID} --destination-path ${year}/${caseID}
    else
	rsync -avz $wkdir/cases/${caseID} answerbe@198.215.54.71:/swnas/cases/answer_devel
    fi
fi
