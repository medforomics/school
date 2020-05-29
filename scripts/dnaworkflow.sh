#!/bin/bash
#SBATCH --job-name dnaworkflow
#SBATCH -N 1
#SBATCH -t 14-0:0:00
#SBATCH -o dnaworkflow.out
#SBATCH -e dnaworkflow.err


usage() {
    echo "-h Help documentation for gatkrunner.sh"
    echo "-p  --projectID"
    echo "-r  --Reference Genome: GRCh38 or GRCm38"
    echo "Example: bash init_workflows.sh -p prefix -r /path/GRCh38"
    exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :r:e:a:p:b:d:m:w:g:c:h opt
do
    case $opt in
        r) index_path=$OPTARG;;
	e) codedir=$OPTARG;;
	a) prodir=$OPTARG;;
        p) prjid=$OPTARG;;
	b) capturedir=$OPTARG;;
	d) mdup=$OPTARG;;
	m) mutectpon=$OPTARG;;
	g) gittag=$OPTARG;;
	h) usage;;
    esac
done

shift $(($OPTIND -1))
baseDir="`dirname \"$0\"`"

if [[ -z $index_path ]]
then
    index_path="/project/shared/bicf_workflow_ref/human/grch38_cloud/dnaref"
fi

outdir="$prodir/fastq"
outnf="$prodir/analysis"
workdir="$prodir/work"

source /etc/profile.d/modules.sh
module load nextflow/20.01.0

pon_opt=''
if [[ -f $mutectpon ]]
then
    pon_opt="--pon $mutectpon"
fi
capture="${capturedir}/targetpanel.bed"

nextflow -C ${codedir}/nextflow.config run -w $workdir ${codedir}/alignment.nf --design design.txt --capturedir $capturedir --capture $capture --input $outdir --output $outnf --markdups $mdup --version $gittag > nextflow_alignment.log

awk '{print "mv '${outnf}'/"$1".* '${outnf}'/"$1"_* '${outnf}'/"$2"/"$1}' design.txt | grep -v SampleID |sh

ln -s ${outnf}/*/*/*.bam $outnf

numsamps=`wc -l design_tumor_only.txt|cut -f 1 -d ' '`
thresh=1
if [[ -f design.txt && "$numsamps" -gt "$thresh" ]]
then
    mkdir tonly
    cd tonly
    nextflow -C ${codedir}/nextflow.config run -w $workdir ${codedir}/tumoronly.nf --design ../design_tumor_only.txt --projectid ${prjid} --input $outnf --output $outnf $pon_opt > nextflow_tumoronly.log &
    cd ..
fi
sleep 20
if [[ -f design_tumor_normal.txt ]]
then
    mkdir somatic
    cd somatic
    nextflow -C ${codedir}/nextflow.config run -w $workdir ${codedir}/somatic.nf --design ../design_tumor_normal.txt --projectid ${prjid} --input $outnf --output $outnf $pon_opt > nextflow_somatic.log &
    cd ..
fi
wait
if [[ -f "${outnf}/GM12878*/GM12878*_${prjid}.dna.vcf.gz" ]]
then
    bamfile=`grep GM12878 design_tumor_only.txt |awk '{print ${outnf}/GM12878/"$1"/"$4}'`
    sampid=`grep GM12878 design_tumor_only.txt |awk '{print $1}'`
    vcfid=`grep GM12878 design_tumor_only.txt |awk '{print $3}'`
    cd ${outnf}/GM12878
    bash $codeir\/scripts/snsp.sh -p $sampid -r $capturedir -t $capture -b $bamfile -v ${vcfid}.dna.vcf.gz
fi
