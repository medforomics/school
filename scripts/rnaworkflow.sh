#!/bin/bash
#SBATCH --job-name rnaworkflow
#SBATCH -N 1
#SBATCH -t 14-0:0:00
#SBATCH -o rnaworkflow.out
#SBATCH -e rnaworkflow.err


usage() {
    echo "-h Help documentation for gatkrunner.sh"
    echo "-p  --projectID"
    echo "-r  --Reference Genome: GRCh38 or GRCm38"
    echo "Example: bash init_workflows.sh -p prefix -r /path/GRCh38"
    exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :r:e:a:p:w:l:c:h opt
do
    case $opt in
        r) index_path=$OPTARG;;
	c) ctbam=$OPTARG;;
	e) codedir=$OPTARG;;
	a) analdir=$OPTARG;;
	l) genelist=$OPTARG;;
        p) prjid=$OPTARG;;
	h) usage;;
    esac
done

shift $(($OPTIND -1))
baseDir="`dirname \"$0\"`"

bamct=''
if [[ -a $ctbam ]]
then
    bamct='--bamct skip'
fi

outdir="$analdir/fastq"
outnf="$analdir/analysis"
workdir="$analdir/work"

module load nextflow/0.31.0
nextflow -C ${codedir}/nextflow.config run -w $workdir ${codedir}/rnaseq.nf --design design.txt --input $outdir --output $outnf $bamct --markdups skip > nextflow.log
ln -fs ${outnf}/*/*/*.bam $outnf
nextflow -C ${codedir}/nextflow.config run -w $workdir ${codedir}/tumoronly.nf --design design_tumor_only.txt --genome $index_path --nuctype rna --callsvs skip --projectid ${prjid} --input $outnf --output $outnf > nextflow_tumoronly.log

if [[ -a $genelist ]]
then
    for i in ${outnf}/*/*/*.fpkm.txt; do
	${codedir}/scripts/fpkm_subset_panel.pl -f $i -g $genelist
    done
fi
