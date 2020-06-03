#!/bin/bash
#init_workflows.sh

set -e

usage() {
  echo "-h Help documentation for gatkrunner.sh"
  echo "-p  --projectID"
  echo "-r  --Reference Genome: GRCh38 or GRCm38"
  echo "-b  --baseDir script executable directory"
  echo "-c  --processing dir"
  echo "-x  --checkXML"
  echo "Example: bash init_workflows.sh -p prefix -r /path/GRCh38 -b baseDir -c processingDir"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :r:p:t:b:s:c:x:h opt
do
    case $opt in
        r) index_path=$OPTARG;;
        p) seqrunid=$OPTARG;;
        c) prodir=$OPTARG;;
	b) baseDir=$OPTARG;;
	x) skipxml=1;;
	t) testing=1;;
	s) splitsamp=1;;
	h) usage;;
    esac
done

shift $(($OPTIND -1))


if [[ -z $baseDir ]]
then
    baseDir="/project/PHG/PHG_Clinical/devel/clinseq_workflows"
fi

cd $baseDir
gittag=`git describe --abbrev=0 --tags`

if [[ -z $index_path ]]
then
    index_path="/project/shared/bicf_workflow_ref/human/grch38_cloud/dnaref"
fi

dna_ref_path=$index_path
rna_ref_path="/project/shared/bicf_workflow_ref/human/grch38_cloud/rnaref"

panelsdir="/project/shared/bicf_workflow_ref/human/grch38_cloud/panels"
illumina="/project/PHG/PHG_Illumina/BioCenter"
fqout="/project/PHG/PHG_Clinical/illumina"

seqdatadir="${fqout}/${seqrunid}"
oriss="${fqout}/sample_sheets/$seqrunid.csv"
newss="${seqdatadir}/$seqrunid.csv"

if [[ -z $prodir ]]
then
    prodir="/project/PHG/PHG_Clinical/processing";
fi
procbase="$prodir/$seqrunid";
outnf="$procbase/analysis"
workdir="$procbase/work"
year="20${seqrunid:0:2}"
mkdir -p /archive/PHG/PHG_Clinical/toarchive/backups/${year}

source /etc/profile.d/modules.sh
mdup='fgbio_umi'
perl ${baseDir}/scripts/create_samplesheet_designfiles.pl -i $oriss -o $newss -d ${prodir}/${seqrunid} -p ${seqrunid} -f ${fqout} -n ${outnf} -t ${panelsdir} -u ${seqdatadir}/$seqrunid.noumi.csv

lastline=`tail -n 1 ${seqdatadir}/$seqrunid.noumi.csv |grep -v Sample_ID`

if [[ -z $testing ]]
then
    module load bcl2fastq/2.19.1
    mkdir -p ${fqout}/${seqrunid}
    ln -s ${illumina}/${seqrunid}/* $fqout/${seqrunid}
    mkdir -p /project/PHG/PHG_BarTender/bioinformatics/demultiplexing/$seqrunid
    if [[ ${#lastline} > 0 ]]
    then
	mkdir -p /project/PHG/PHG_BarTender/bioinformatics/demultiplexing/$seqrunid/noumi
	bcl2fastq --barcode-mismatches 0 -o ${seqdatadir} --no-lane-splitting --runfolder-dir ${seqdatadir} --sample-sheet ${seqdatadir}/$seqrunid.noumi.csv --use-bases-mask Y76,I6N8,Y76 &> ${seqdatadir}/bcl2fastq_noumi_${seqrunid}.log
	rsync -avz ${seqdatadir}/Reports ${seqdatadir}/Stats answerbe@198.215.54.71:/swnas/qc_nuclia/demultiplexing/$seqrunid/noumi
	mkdir -p ${seqdatadir}/noumi
	mv ${seqdatadir}/Reports ${seqdatadir}/Stats ${seqdatadir}/noumi
    fi
    mv ${seqdatadir}/RunInfo.xml ${seqdatadir}/RunInfo.xml.ori
    perl  ${baseDir}/scripts/fix_runinfo_xml.pl $seqdatadir
    
    bcl2fastq --barcode-mismatches 0 -o ${seqdatadir} --no-lane-splitting --runfolder-dir ${seqdatadir} --sample-sheet ${newss} &> ${seqdatadir}/bcl2fastq_${seqrunid}.log
    rsync -avz ${seqdatadir}/Reports ${seqdatadir}/Stats answerbe@198.215.54.71:/swnas/qc_nuclia/demultiplexing/$seqrunid
    rsync -avz --prune-empty-dirs --include "*/" --include="*.fastq.gz" --include="*.csv" --exclude="*" $seqdatadir /archive/PHG/PHG_Clinical/toarchive/seqruns/$year
    source ${baseDir}/../azure_credentials
    az storage blob upload-batch -d seqruns -s /archive/PHG/PHG_Clinical/toarchive/seqruns/${year}/${seqrunid} --destination-path ${year}/${seqrunid}
fi

bash ${prodir}/${seqrunid}/lnfq.sh 
cd ${prodir}/${seqrunid}

panelbed=(["panel1385"]="UTSW_V2_pancancer" ["panel1385v2"]="UTSW_V3_pancancer" ["medexomeplus"]="UTSW_medexomeplus" ["heme"]="UTSW_V4_heme" ["pancancer"]="UTSW_V4_pancancer" ["idtrnaseq"]="UTSW_V4_rnaseq" ["solid"]="UTSW_V4_solid" ["tspcrnaseq"]="TSPC" )

while read i; do
    caseid=$i
    cd ${prodir}/${seqrunid}/analysis/$caseid
    if [[ -z $skipxml ]] && [[ ! -f $i.xml ]] && [[ -n $clarityuser ]]
    then
	python ${baseDir}/IntellispaceDemographics/gatherdemographics.py -i $i -u $clarityuser -p $claritypasswd -o ${i}.xml
	missing=`grep '><\|D64.9\|N\/A' ${i}.xml`
	if [[ -a $missing ]]
	then
	    SUBJECT="SECUREMAIL: Case Missing Data"
	    TO="ngsclialab@UTSouthwestern.edu,Hui.Zheng@UTSouthwestern.edu,Yan.Xu@UTSouthwestern.edu"
	    email=$baseDir"/scripts/bioinformatics_email.txt"
	    cat $email | sed 's/Unspecified/'$i'/' | /bin/mail -s "$SUBJECT" "$TO"
	fi
    fi
done <subjects.txt

if [[ -z $splitsamp ]]
then
    for i in */design.txt; do
	dtype=`dirname $i`
	cd ${prodir}/${seqrunid}/${dtype}
	bash lnfq.sh  ###remove symlinks to analysis
	if [[ $dtype == 'wholernaseq' ]]
	then
	    cp ${baseDir}/scripts/rnaworkflow.sh ${procbase}
	    bash ${procbase}/rnaworkflow.sh -r $rna_ref_path -e $baseDir -a $procbase -p $seqrunid -c 1 -g $gittag &> log.txt &
	elif [[ $dtype =~ 'rnaseq' ]]
	then
	    genelist="${panelsdir}/${panelbed[${dtype}]}/genelist.txt"
	    cp ${baseDir}/scripts/rnaworkflow.sh ${procbase}
	    bash ${procbase}/rnaworkflow.sh -r $rna_ref_path -e $baseDir -a $procbase -p $seqrunid -l $genelist -g $gittag &> log.txt &
	else
	    pon_opt=''
	    if [[ -f "${panelsdir}/${panelbed[${dtype}]}/mutect2.pon.vcf.gz" ]]
	    then
		pon_opt="-m ${panelsdir}/${panelbed[${dtype}]}/mutect2.pon.vcf.gz"
	    fi
	    capture="${panelsdir}/${panelbed[${dtype}]}"
	    cp ${baseDir}/scripts/dnaworkflow.sh ${procbase}
	    bash ${procbase}/dnaworkflow.sh -r $dna_ref_path -e $baseDir -b $capture -a $procbase -p $seqrunid -d $mdup $pon_opt -g $gittag &> log.txt &
	fi
    done
    wait
    cd $outnf
    cd $prodir\/$seqrunid
    if [[ -z $testing ]]
    then
	rsync -rlptgoD --exclude="*fastq.gz*" --exclude "*work*" --exclude="*bam*" ${prodir}/${seqrunid} answerbe@198.215.54.71:/swnas/qc_nuclia/seqanalysis
	for i in *.properties; do
	    curl "http://nuclia.biohpc.swmed.edu:8080/NuCLIAVault/addPipelineResultsWithProp?token=${nucliatoken}&propFilePath=/swnas/qc_nuclia/seqanalysis/$seqrunid/${i}"
	done
	rsync -avz --prune-empty-dirs --include "*/" --include="*.fastq.gz" --include="*.csv" --exclude="*" $seqdatadir /archive/PHG/PHG_Clinical/toarchive/seqruns/$year
	rsync -avz --no-links --exclude="*fastq.gz*" ${prodir}/${seqrunid}/analysis/* /archive/PHG/PHG_Clinical/cases/
	source ${baseDir}/../azure_credentials
	az storage blob upload-batch -d seqruns -s /archive/PHG/PHG_Clinical/toarchive/seqruns/${year}/${seqrunid} --destination-path ${year}/${seqrunid}
    fi
else
    while read i; do
	for j in */vars.sh; do #$dtype/vars.sh
	    dtype=`echo $j |cut -f 1 -d '/'`
	    cd ${prodir}/${seqrunid}/analysis/$caseid/$dtype
	    cat $j > run_wkflow.sh
	    echo "module load nextflow/20.01.0" >> run_wkflow.sh
	    if [[ $dtype == 'wholernaseq' ]]
	    then
		echo "nextflow -C ${baseDir}/nextflow.config run -w $workdir ${baseDir}/rna.nf --input \$input --output \$output --caseid $caseid --seqrunid $seqrunid --tumorid \$tumorid --bamct skip --version $gittag --genome $rna_ref_path -resume" >> run_wkflow.sh
	    elif [[ $dtype =~ 'rna' ]]
	    then
		genelist="${panelsdir}/${panelbed[${dtype}]}/genelist.txt"
		echo "nextflow -C ${baseDir}/nextflow.config run -w $workdir ${baseDir}/rna.nf --input \$input --output \$output --caseid $caseid --seqrunid $seqrunid --tumorid \$tumorid --version $gittag  --genome $rna_ref_path --glist $genelist -resume" >> run_wkflow.sh
	    else
		normalopt=''
		if [[ `grep _N_ samples.txt` ]]
		then
		    normalopt="--normalid \$normalid"
		fi
		echo "nextflow -C ${baseDir}/nextflow.config run -w $workdir ${baseDir}/dna.nf --input \$input --output \$output --caseid $caseid --seqrunid $seqrunid --tumorid \$tumorid --pon \$mutectpon --capture \$capturebed --capturedir \$capturedir --version $gittag --genome $dna_ref_path $normalopt -resume" >> run_wkflow.sh
		if [[ ${caseid} =~ "GM12878" ]]
		then
		    echo "cd ${prodir}/${seqrunid}/analysis/$caseid" >> run_wkflow.sh
		    echo "bash ${baseDir}/scripts/snsp.sh -p \$sampid -t \$capturebed -b \$bamfile -v \$vcf" >> run_wkflow.sh
		fi
	    fi
	    if [[ -z $testing ]]
	    then
		echo "rsync -rlptgoD --exclude=\"*fastq.gz*\" --exclude \"*work*\" --exclude=\"*bam*\" ${prodir}/${seqrunid} answerbe@198.215.54.71:/swnas/qc_nuclia/seqanalysis"  >> run_wkflow.sh
		echo "rsync -avz --no-links --exclude=\"*fastq.gz*\" ${prodir}/${seqrunid}/analysis/${caseid} /archive/PHG/PHG_Clinical/cases/"  >> run_wkflow.sh
		while read j; do
		    echo "curl \"http://nuclia.biohpc.swmed.edu:8080/NuCLIAVault/addPipelineResultsWithProp?token=${nucliatoken}&propFilePath=/swnas/qc_nuclia/seqanalysis/${seqrunid}/${j}.properties\"" >> run_wkflow.sh
		done < samples.txt
	    fi
	    sbatch -p 32GB,super run_wkflow.sh
	done
    done <subjects.txt
fi
