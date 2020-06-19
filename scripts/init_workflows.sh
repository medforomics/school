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
while getopts :r:p:t:b:c:x:h opt
do
    case $opt in
        r) index_path=$OPTARG;;
        p) seqrunid=$OPTARG;;
        c) prodir=$OPTARG;;
	b) baseDir=$OPTARG;;
	x) skipxml=1;;
	t) testing=1;;
	h) usage;;
    esac
done

shift $(($OPTIND -1))

echo "*****Setting Variables******"
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
echo "*****DONE Setting Variables******"


echo "*****Creating Samplesheets******"
source /etc/profile.d/modules.sh
mdup='fgbio_umi'
mkdir -p ${fqout}/${seqrunid}
echo "${baseDir}/scripts/create_samplesheet_designfiles.pl -i $oriss -o $newss -d ${prodir}/${seqrunid} -p ${seqrunid} -f ${fqout} -n ${outnf} -t ${panelsdir} -u ${seqdatadir}/$seqrunid.noumi.csv"
perl ${baseDir}/scripts/create_samplesheet_designfiles.pl -i $oriss -o $newss -d ${prodir}/${seqrunid} -p ${seqrunid} -f ${fqout} -n ${outnf} -t ${panelsdir} -u ${seqdatadir}/$seqrunid.noumi.csv
echo "*****Done Creating Samplesheets******"

echo "*****Starting Demultiplexing******"
lastline=`tail -n 1 ${seqdatadir}/$seqrunid.noumi.csv |cut -f 1 -d ','`
echo $lastline

if [[ -z $testing ]]
then
    module load bcl2fastq/2.19.1
    ln -s ${illumina}/${seqrunid}/* $fqout/${seqrunid}
    ssh answerbe@198.215.54.71 "mkdir -p /swnas/qc_nuclia/demultiplexing/$prjid"
    if [[ $lastline != SampleID ]]
    then
	ssh answerbe@198.215.54.71 "mkdir -p /swnas/qc_nuclia/demultiplexing/$prjid/noumi"
	bcl2fastq --barcode-mismatches 0 -o ${seqdatadir} --no-lane-splitting --runfolder-dir ${seqdatadir} --sample-sheet ${seqdatadir}/$prjid.noumi.csv --use-bases-mask Y76,I6N
	8,Y76 &> ${seqdatadir}/bcl2fastq_noumi_${prjid}.log
	if [[ ! -s "${seqdatadir}/Undetermined_S0_R1_001.fastq.gz" ]] || [[ ! -s "${seqdatadir}/Undetermined_S0_R2_001.fastq.gz" ]]
	then
	    exit
	fi
	rsync -avz ${seqdatadir}/Reports ${seqdatadir}/Stats answerbe@198.215.54.71:/swnas/qc_nuclia/demultiplexing/$prjid/noumi
	mkdir -p ${seqdatadir}/noumi
	mv ${seqdatadir}/Reports ${seqdatadir}/Stats ${seqdatadir}/noumi
    fi
    mv ${seqdatadir}/RunInfo.xml ${seqdatadir}/RunInfo.xml.ori
    perl  ${baseDir}/scripts/fix_runinfo_xml.pl $seqdatadir
    bcl2fastq --barcode-mismatches 0 -o ${seqdatadir} --no-lane-splitting --runfolder-dir ${seqdatadir} --sample-sheet ${newss} &> ${seqdatadir}/bcl2fastq_${prjid}.log
    if [[ ! -s "${seqdatadir}/Undetermined_S0_R1_001.fastq.gz" ]] || [[ ! -s "${seqdatadir}/Undetermined_S0_R2_001.fastq.gz" ]]
    then
	exit
    fi
    rsync -avz ${seqdatadir}/Reports ${seqdatadir}/Stats answerbe@198.215.54.71:/swnas/qc_nuclia/demultiplexing/$prjid/
    rsync -avz --prune-empty-dirs --include "*/" --include="*.fastq.gz" --include="*.csv" --exclude="*" $seqdatadir /archive/PHG/PHG_Clinical/toarchive/seqruns/$year
    source ${baseDir}/../azure_credentials
    az storage blob upload-batch -d seqruns -s /archive/PHG/PHG_Clinical/toarchive/seqruns/${year}/${seqrunid} --destination-path ${year}/${seqrunid}
fi

echo "*****Done Demultiplexing******"

cd ${prodir}/${seqrunid}
bash lnfq.sh

panelbed=(["panel1385"]="UTSW_V2_pancancer" ["panel1385v2"]="UTSW_V3_pancancer" ["medexomeplus"]="UTSW_medexomeplus" ["heme"]="UTSW_V4_heme" ["pancancer"]="UTSW_V4_pancancer" ["idtrnaseq"]="UTSW_V4_rnaseq" ["solid"]="UTSW_V4_solid" ["tspcrnaseq"]="TSPC" )

echo "*****Checking Case Information******"

while read i; do
    caseid=$i
#    cd ${prodir}/${seqrunid}/analysis/$caseid
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

echo "*****Creating run_workflow bash scripts******"

for i in */design.txt; do
    dtype=`dirname $i`
    cd ${prodir}/${seqrunid}/${dtype}
    cat vars.sh > run_wkflow.sh
    echo "module load nextflow/20.01.0" >> run_wkflow.sh
    if [[ $dtype == 'wholernaseq' ]]
    then
	echo "nextflow -C ${baseDir}/nextflow.config run -w $workdir ${baseDir}/rna.nf --input \$input --output \$output --seqrunid $seqrunid --design.txt --bamct skip --version $gittag --genome $rna_ref_path -resume" >> run_wkflow.sh
    elif [[ $dtype =~ 'rna' ]]
    then
	genelist="${panelsdir}/${panelbed[${dtype}]}/genelist.txt"
	echo "nextflow -C ${baseDir}/nextflow.config run -w $workdir ${baseDir}/rna.nf --input \$input --output \$output --design design.txt --seqrunid $seqrunid --version $gittag  --genome $rna_ref_path --glist $genelist -resume" >> run_wkflow.sh
    else
	echo "nextflow -C ${baseDir}/nextflow.config run -w $workdir ${baseDir}/dna.nf --input \$input --output \$output --design design.txt --seqrunid $seqrunid --pon \$mutectpon --capture \$capturebed --capturedir \$capturedir --version $gittag --genome $dna_ref_path -resume" >> run_wkflow.sh
	if [[ `grep GM12878 design.txt` ]]
	then
	    caseid=`grep GM12878 design.txt |cut -f 2`
	    echo "cd ${prodir}/${seqrunid}/analysis/$caseid" >> run_wkflow.sh
	    echo "bash ${baseDir}/scripts/snsp.sh -p \$sampid -t \$capturebed -b \$bamfile -v \$vcf" >> run_wkflow.sh
	fi
    fi
    if [[ -z $testing ]]
    then
	echo "rsync -rlptgoD --exclude=\"*fastq.gz*\" --exclude \"*work*\" --exclude=\"*bam*\" ${prodir}/${seqrunid} answerbe@198.215.54.71:/swnas/qc_nuclia/seqanalysis"  >> run_wkflow.sh
	echo "rsync -avz --no-links --exclude=\"*fastq.gz*\" ${prodir}/${seqrunid}/analysis/${caseid} /archive/PHG/PHG_Clinical/cases/"  >> run_wkflow.sh
	echo "cd ${prodir}/${seqrunid}" >> run_wkflow.sh
	for j in *properties ; do
	    echo "curl \"http://nuclia.biohpc.swmed.edu:8080/NuCLIAVault/addPipelineResultsWithProp?token=\${nucliatoken}&propFilePath=/swnas/qc_nuclia/seqanalysis/${seqrunid}/${j}\"" >> run_wkflow.sh
	done
    fi
    sbatch -p 32GB,super run_wkflow.sh
done
