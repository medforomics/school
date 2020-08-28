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
    baseDir="/project/PHG/PHG_Clinical/devel/school"
fi
pendir=/project/PHG/PHG_Clinical/cloud/pending
cd $baseDir

illumina="/project/PHG/PHG_Illumina/BioCenter"
fqout="/project/PHG/PHG_Clinical/illumina"
seqdatadir="${fqout}/${seqrunid}"
oriss="${fqout}/sample_sheets/$seqrunid.csv"
newss="${seqdatadir}/$seqrunid.csv"

if [[ -z $prodir ]]
then
    prodir="/project/PHG/PHG_Clinical/devel/processing";
fi
procbase="$prodir/$seqrunid";
outnf="$procbase/analysis"
year="20${seqrunid:0:2}"
echo "*****DONE Setting Variables******"


echo "*****Creating Samplesheets******"
source /etc/profile.d/modules.sh
module load perl/5.28.0
mdup='fgbio_umi'
mkdir -p ${prodir}/${seqrunid}
mkdir -p ${fqout}/${seqrunid}
perl ${baseDir}/scripts/update_samplesheet.pl -i $oriss -o $newss
perl ${baseDir}/scripts/create_redmine_issue.pl $newss
echo "*****Done Creating Samplesheets******"

echo "*****Starting Demultiplexing******"

if [[ -z $testing ]]
then
    mkdir -p /archive/PHG/PHG_Clinical/toarchive/backups/${year}
    module load bcl2fastq/2.19.1
    ln -s ${illumina}/${seqrunid}/* $fqout/${seqrunid}
    ssh answerbe@198.215.54.71 "mkdir -p /swnas/qc_nuclia/demultiplexing/$prjid"
    mv ${seqdatadir}/RunInfo.xml ${seqdatadir}/RunInfo.xml.ori
    perl  ${baseDir}/scripts/fix_runinfo_xml.pl $seqdatadir
    bcl2fastq --barcode-mismatches 0 -o ${seqdatadir} --no-lane-splitting --runfolder-dir ${seqdatadir} --sample-sheet ${newss} &> ${seqdatadir}/bcl2fastq_${prjid}.log
    if [[ ! -s "${seqdatadir}/Undetermined_S0_R1_001.fastq.gz" ]] || [[ ! -s "${seqdatadir}/Undetermined_S0_R2_001.fastq.gz" ]]
    then
	exit
    fi
    rsync -avz ${seqdatadir}/Reports ${seqdatadir}/Stats answerbe@198.215.54.71:/swnas/qc_nuclia/demultiplexing/$prjid/
    rsync -avz --prune-empty-dirs --include "*/" --include="*.fastq.gz" --include="*.csv" --exclude="*" $seqdatadir /archive/PHG/PHG_Clinical/toarchive/seqruns/$year
    source ${baseDir}/azure_credentials
    az storage blob upload-batch -d seqruns -s /archive/PHG/PHG_Clinical/toarchive/seqruns/${year}/${seqrunid} --destination-path ${year}/${seqrunid}
fi

echo "*****Done Demultiplexing******"

cd ${prodir}/${seqrunid}
mkdir fastq
find $seqdatadir -name "*.fastq.gz" -exec ln -s {} fastq \;
find $seqdatadir -name "*.fastq.gz" -print |perl -pe 's/\//\t/g' |cut -f 6- |sort |grep -v Undetermined > sequence.files.txt
cut -f 2 sequence.files.txt | sort -u > subjects.txt
perl ${baseDir}/scripts/create_designfiles.pl -i $newss -o sequence.files.txt

echo "*****Submitting jobs******"
for df in *.design.txt
do
    wkflow="${df%.design.txt}"
    echo "bash ${appletdir}/workflows/run_workflow.sh -d $df -w $wkflow -i fastq -p $seqrunid -o $pendir"
done

echo "*****Checking Case Information******"

while read i; do
    caseid=$i
    if [[ -z $skipxml ]] && [[ ! -f $i.xml ]] && [[ -n $clarityuser ]]
    then
	python ${baseDir}/IntellispaceDemographics/gatherdemographics.py -i $i -u $clarityuser -p $claritypasswd -o ${i}.xml
	missing=`grep '><\|D64.9\|N\/A' ${i}.xml`
	if [[ -a $missing ]]
	then
	    SUBJECT="SECUREMAIL: Case Missing Data"
	    TO="ngsclialab@UTSouthwestern.edu,Hui.Zheng@UTSouthwestern.edu,Yan.Xu@UTSouthwestern.edu"x
	    email=$baseDir"/scripts/bioinformatics_email.txt"
	    cat $email | sed 's/Unspecified/'$i'/' | /bin/mail -s "$SUBJECT" "$TO"
	fi
    fi
done <subjects.txt


