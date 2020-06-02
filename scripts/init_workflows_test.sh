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
while getopts :r:p:b:c:x:h opt
do
    case $opt in
        r) index_path=$OPTARG;;
        p) prjid=$OPTARG;;
        c) prodir=$OPTARG;;
	b) baseDir=$OPTARG;;
	x) skipxml=1;;
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

seqdatadir="${fqout}/${prjid}"
oriss="${fqout}/sample_sheets/$prjid.csv"
newss="${seqdatadir}/$prjid.csv"

if [[ -z $prodir ]]
then
    prodir="/project/PHG/PHG_Clinical/processing";
fi
procbase="$prodir/$prjid";
outdir="$procbase/fastq"
outnf="$procbase/analysis"
workdir="$procbase/work"

mkdir -p ${fqout}/${prjid}
ln -s ${illumina}/${prjid}/* $fqout/${prjid}

declare -A panelbed
panelbed=(["panel1385"]="UTSW_V2_pancancer" ["panel1385v2"]="UTSW_V3_pancancer" ["medexomeplus"]="UTSW_medexomeplus" ["heme"]="UTSW_V4_heme" ["pancancer"]="UTSW_V4_pancancer" ["idtrnaseq"]="UTSW_V4_rnaseq" ["solid"]="UTSW_V4_solid")

source /etc/profile.d/modules.sh
module load bcl2fastq/2.19.1

#mkdir -p /project/PHG/PHG_BarTender/bioinformatics/demultiplexing/$prjid
perl ${baseDir}/scripts/create_samplesheet_designfiles.pl -i $oriss -o $newss -d ${prodir}/${prjid} -p ${prjid} -f ${outdir} -n ${outnf} -u ${seqdatadir}/$prjid.noumi.csv
lastline=`tail -n 1 ${seqdatadir}/$prjid.noumi.csv |grep -v Sample_ID`

if [[ ${#lastline} > 0 ]]
then
#    mkdir -p /project/PHG/PHG_BarTender/bioinformatics/demultiplexing/$prjid/noumi
#    bcl2fastq --barcode-mismatches 0 -o ${seqdatadir} --no-lane-splitting --runfolder-dir ${seqdatadir} --sample-sheet ${seqdatadir}/$prjid.noumi.csv --use-bases-mask Y76,I6N8,Y76 &> ${seqdatadir}/bcl2fastq_noumi_${prjid}.log
#    rsync -avz ${seqdatadir}/Reports ${seqdatadir}/Stats answerbe@198.215.54.71:/swnas/qc_nuclia/demultiplexing/$prjid/noumi
#    mkdir -p ${seqdatadir}/noumi
#    mv ${seqdatadir}/Reports ${seqdatadir}/Stats ${seqdatadir}/noumi
fi
#mv ${seqdatadir}/RunInfo.xml ${seqdatadir}/RunInfo.xml.ori
#perl  ${baseDir}/scripts/fix_runinfo_xml.pl $seqdatadir
mdup='fgbio_umi'

#bcl2fastq --barcode-mismatches 0 -o ${seqdatadir} --no-lane-splitting --runfolder-dir ${seqdatadir} --sample-sheet ${newss} &> ${seqdatadir}/bcl2fastq_${prjid}.log
#rsync -avz ${seqdatadir}/Reports ${seqdatadir}/Stats answerbe@198.215.54.71:/swnas/qc_nuclia/demultiplexing/$prjid
cd ${prodir}/${prjid}

for i in */design.txt; do
    dtype=`dirname $i`
    cd ${prodir}/${prjid}/${dtype}
    subjs=`cut -f 2 design.txt |grep -v FamilyID |uniq`
    for i in $subjs; do
	mkdir -p $outnf/$i
	mkdir -p $outnf/$i/fastq
	if [[ -z $skipxml ]]
	then
	    if [[ ! -f $i.xml ]]
	    then
		python ${baseDir}/IntellispaceDemographics/gatherdemographics.py -i $i -u phg_workflow -p ${password} -o ${i}.xml
		missing=`grep '><\|D64.9\|N\/A' ${i}.xml`
		if [[ -n $missing ]]
		then
		    SUBJECT="SECUREMAIL: Case Missing Data"
		    TO="erika.villa@utsouthwestern.edu,Hui.Zheng@UTSouthwestern.edu,Yan.Xu@UTSouthwestern.edu"
		    email=$baseDir"/scripts/bioinformatics_email.txt"
		    cat $email | sed 's/Unspecified/'$i'/' | /bin/mail -s "$SUBJECT" "$TO"
		fi
	    fi
	fi
    done
    bash lnfq.sh 
    if [[ $dtype == 'wholernaseq' ]]
    then
	cp ${baseDir}/scripts/rnaworkflow.sh ${procbase}
	bash ${procbase}/rnaworkflow.sh -r $rna_ref_path -e $baseDir -a $procbase -p $prjid -c 1 -g $gittag &> log.txt &
    elif [[ $dtype =~ 'rnaseq' ]]
    then
	genelist="${panelsdir}/${panelbed[${dtype}]}/genelist.txt"
	cp ${baseDir}/scripts/rnaworkflow.sh ${procbase}
	bash ${procbase}/rnaworkflow.sh -r $rna_ref_path -e $baseDir -a $procbase -p $prjid -l $genelist -g $gittag &> log.txt &
    else
	pon_opt=''
	if [[ -f "${panelsdir}/${panelbed[${dtype}]}/mutect2.pon.vcf.gz" ]]
	then
	    pon_opt="-m ${panelsdir}/${panelbed[${dtype}]}/mutect2.pon.vcf.gz"
	fi
	capture="${panelsdir}/${panelbed[${dtype}]}"
	cp ${baseDir}/scripts/dnaworkflow.sh ${procbase}
	bash ${procbase}/dnaworkflow.sh -r $dna_ref_path -e $baseDir -b $capture -a $procbase -p $prjid -d $mdup $pon_opt -g $gittag &> log.txt &
    fi
done
wait
cd $outnf

# cd $prodir\/$prjid
# rsync -rlptgoD --exclude="*fastq.gz*" --exclude "*work*" --exclude="*bam*" ${prodir}/${prjid} answerbe@198.215.54.71:/swnas/qc_nuclia/seqanalysis
# #perl ${baseDir}/scripts/create_properties_run.pl -p $prjid -d $prodir
# #rsync -rlptgoD --exclude="*fastq.gz*" --exclude "*work*" --exclude="*bam*" ${prodir}/${prjid} answerbe@198.215.54.71:/swnas/qc_nuclia/seqanalysis

# #for i in *.properties; do
# #    curl "http://nuclia.biohpc.swmed.edu:8080/NuCLIAVault/addPipelineResultsWithProp?token=${nucliatoken}&propFilePath=/swnas/qc_nuclia/seqanalysis/$prjid/${i}"
# done

# year="20${prjid:0:2}"
# mkdir -p /archive/PHG/PHG_Clinical/toarchive/backups/${year}

# rsync -avz --prune-empty-dirs --include "*/" --include="*.fastq.gz" --include="*.csv" --exclude="*" $seqdatadir /archive/PHG/PHG_Clinical/toarchive/seqruns/$year
# rsync -avz --no-links --exclude="*fastq.gz*" ${prodir}/${prjid}/analysis/* /archive/PHG/PHG_Clinical/cases/

# source ${baseDir}/../azure_credentials
# az storage blob upload-batch -d seqruns -s /archive/PHG/PHG_Clinical/toarchive/seqruns/${year}/${prjid} --destination-path ${year}/${prjid}
