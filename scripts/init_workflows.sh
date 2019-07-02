#!/bin/bash
#init_workflows.sh

usage() {
  echo "-h Help documentation for gatkrunner.sh"
  echo "-p  --projectID"
  echo "-r  --Reference Genome: GRCh38 or GRCm38"
  echo "Example: bash init_workflows.sh -p prefix -r /path/GRCh38"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :r:p:c::xh opt
do
    case $opt in
        r) index_path=$OPTARG;;
        p) prjid=$OPTARG;;
        c) prodir=$OPTARG;;
	x) checkxml=1;;
	h) usage;;
    esac
done

shift $(($OPTIND -1))
baseDir="`dirname \"$0\"`"

fqout="/project/PHG/PHG_Clinical/illumina"
clinref="/project/shared/bicf_workflow_ref/human/GRCh38/clinseq_prj"
hisat_index_path='/project/shared/bicf_workflow_ref/human/GRCh38/hisat_index'
illumina="/project/PHG/PHG_Illumina/BioCenter"
if [[ -z $prodir ]]
then
prodir="/project/PHG/PHG_Clinical/processing";
fi

seqdatadir="${fqout}/${prjid}"
oriss="${fqout}/sample_sheets/$prjid.csv"
newss="${seqdatadir}/$prjid.csv"

mkdir ${fqout}/${prjid}
ln -s ${illumina}/${prjid}/* $fqout/${prjid}

umi=`grep "<Read Number=\"2\" NumCycles=\"14\" IsIndexedRead=\"Y\" />" ${illumina}/${prjid}/RunInfo.xml`

outdir="$prodir/$prjid/fastq"
outnf="$prodir/$prjid/analysis"
workdir="$prodir/$prjid/work"

if [[ ! -d "$prodir/$prjid" ]]
then 
    mkdir ${prodir}/${prjid}
    mkdir $outdir
    mkdir $outnf
    mkdir $workdir
fi
declare -A panelbed
panelbed=(["panel1385"]="UTSWV2.bed" ["panel1385v2"]="UTSWV2_2.panelplus.bed" ["idthemev1"]="heme_panel_probes.bed" ["idthemev2"]="hemepanelV3.bed" ["idtcellfreev1"]="panelcf73_idt.100plus.bed" ["medexomeplus"]="MedExome_Plus.bed")

if [[ -n $umi ]]
then
    mv ${seqdatadir}/RunInfo.xml ${seqdatadir}/RunInfo.xml.ori
    perl $baseDir/fix_runinfo_xml.pl $seqdatadir
    perl $baseDir/create_samplesheet_designfiles.pl -i $oriss -o $newss -d ${prodir}/${prjid} -p ${prjid} -f ${outdir} -n ${outnf} -u
    mdup='fgbio_umi'
else
    perl $baseDir/create_samplesheet_designfiles.pl -i $oriss -o $newss -d ${prodir}/${prjid} -p ${prjid} -f ${outdir} -n ${outnf}
    mdup='picard'
fi

source /etc/profile.d/modules.sh
module load bcl2fastq/2.19.1 nextflow/0.31.0 vcftools/0.1.14 samtools/gcc/1.8
bcl2fastq --barcode-mismatches 0 -o ${seqdatadir} --no-lane-splitting --runfolder-dir ${seqdatadir} --sample-sheet ${newss} &> ${seqdatadir}/bcl2fastq_${prjid}.log
if [[ ! -d /project/PHG/PHG_BarTender/bioinformatics/demultiplexing/${prjid} ]]
then
   mkdir /project/PHG/PHG_BarTender/bioinformatics/demultiplexing/$prjid\n
fi
rsync -avz ${seqdatadir}/Reports /project/PHG/PHG_BarTender/bioinformatics/demultiplexing/$prjid
rsync -avz ${seqdatadir}/Stats /project/PHG/PHG_BarTender/bioinformatics/demultiplexing/$prjid
cd ${prodir}/${prjid}
codedir="/project/PHG/PHG_Clinical/devel/clinseq_workflows"
for i in */design.txt; do
    dtype=`dirname $i`
    cd ${prodir}/${prjid}/${dtype}
    subjs=`cut -f 2 design.txt |grep -v FamilyID |uniq`
    for i in $subjs; do
	mkdir -p $outnf/$i
	mkdir -p $outnf/$i/fastq
	if [[ -a $checkxml ]]
	then
	    if [[ ! -f $i.xml ]]
	    then
		python $codedir/../IntellispaceDemographics/gatherdemographics.py -i $i -u phg_workflow -p $password -o ${i}.xml
		missing=`grep '><\|D64.9\|N\/A' ${i}.xml`
		if [[ -n $missing ]]
		then
		    SUBJECT="SECUREMAIL: Case Missing Data"
		    TO="erika.villa@utsouthwestern.edu,Hui.Zheng@UTSouthwestern.edu,Yan.Xu@UTSouthwestern.edu"
		    email=$codedir"/bioinformatics_email.txt"
		    cat $email | sed 's/Unspecified/'$i'/' | /bin/mail -s "$SUBJECT" "$TO"
		fi
	    fi
	fi
    done
    bash lnfq.sh
    if [[ $dtype == 'panelrnaseq' ]]
    then
	bash ${codedir}/scripts/rnaworkflow.sh -r $hisat_index_path -e $codedir -a "$prodir/$prjid" -p $prjid &> log.txt &
    elif [[ $dtype == 'wholernaseq' ]]
	then
	    bash ${codedir}/scripts/rnaworkflow.sh -r $hisat_index_path -e $codedir -a "$prodir/$prjid" -p $prjid -c &> log.txt &
    else
	capture="${clinref}/${panelbed[${dtype}]}"
	bash ${codedir}/scripts/dnaworkflow.sh -r $index_path -e $codedir -b $capture -a "$prodir/$prjid" -p $prjid -d $mdup &> log.txt &
    fi
done
wait
cd $outnf

cd $prodir\/$prjid
rsync -rlptgoD --exclude="*fastq.gz*" --exclude "*work*" --exclude="*bam*" ${prodir}/${prjid} /project/PHG/PHG_BarTender/bioinformatics/seqanalysis/
perl ${baseDir}/scripts/create_properties_run.pl -p $prjid -d /project/PHG/PHG_BarTender/bioinformatics/seqanalysis

for i in /project/PHG/PHG_BarTender/bioinformatics/seqanalysis/${prjid}/*.properties; do
    curl "http://nuclia.biohpc.swmed.edu:8080/NuCLIAVault/addPipelineResultsWithProp?token=${nucliatoken}&propFilePath=${i}"
done
monyear="20${prjid:0:4}"
if [[ ! -d "/archive/PHG/PHG_Clinical/toarchive/backups/${monyear}" ]]
then
    mkdir /archive/PHG/PHG_Clinical/toarchive/backups/${monyear}
fi

tar cf /work/archive/PHG/PHG_Clinical/toarchive/backups/${monyear}/${prjid}.tar.gz $seqdatadir
