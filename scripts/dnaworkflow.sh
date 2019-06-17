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
while getopts :r:e:a:p:b:d:w:c:h opt
do
    case $opt in
        r) index_path=$OPTARG;;
	e) codedir=$OPTARG;;
	a) prodir=$OPTARG;;
        p) prjid=$OPTARG;;
	b) capture=$OPTARG;;
	d) mdup=$OPTARG;;
	h) usage;;
    esac
done

shift $(($OPTIND -1))
baseDir="`dirname \"$0\"`"

if [[ -z $index_path ]]
then
    index_path="/project/shared/bicf_workflow_ref/human/GRCh38"
fi
capturedir="${index_path}/clinseq_prj"
   
outdir="$prodir/fastq"
outnf="$prodir/analysis"
workdir="$prodir/work"

module load nextflow/0.31.0

nextflow -C ${codedir}/nextflow.config run -w $workdir ${codedir}/alignment.nf --design design.txt --capture $capture --input $outdir --output $outnf --markdups $mdup > nextflow_alignment.log

mv ${outnf}/*/*/*.bam $outnf
awk '{print "mv ${outnf}/"$1".* ${outnf}/"$1"_* ${outnf}/"$2"/"$1}' design.txt | grep -v SampleID |sh
numsamps=`wc -l design_tumor_only.txt|cut -f 1 -d ' '`
thresh=2
if [[ -f design.txt && "$numsamps" -gt "$thresh" ]]
then
    nextflow -C ${codedir}/nextflow.config run -w $workdir ${codedir}/tumoronly.nf --design design_tumor_only.txt --projectid ${prjid} --input $outnf --targetpanel $capture --output $outnf > nextflow_tumoronly.log &
fi

if [[ -f design_tumor_normal.txt ]]
then
    nextflow -C ${codedir}/nextflow.config run -w $workdir ${codedir}/somatic.nf --design design_tumor_normal.txt --projectid ${prjid} --input $outnf --output $outnf > nextflow_somatic.log &
fi
wait
if [[ -f "${outnf}/GM12878/GM12878_dna_${prjid}.germline.vcf.gz" ]]
then
    cd ${outnf}/GM12878
    bamfile=`grep GM12878 design_tumor_only.txt |awk '{print ${outnf}/GM12878/"$1"/"$4}'`
    bash $codeir\/scripts/snsp.sh -p $sampid -r $capturedir -t $capture -b $bamfile -v ${outnf}/GM12878/GM12878_dna_${prjid}.germline.vcf.gz
fi
