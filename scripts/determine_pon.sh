#!/bin/bash
#determine_pon.sh

usage() {
  echo "-h Help documentation for hisat.sh"
  echo "-r  --Reference Genome: GRCh38 or GRCm38"
  echo "-b  --BAM File"
  echo "-c  --Capture Bedfile"
  echo "Example: bash determine_pon.sh -r /project/shared/bicf_workflow_ref/human/GRCh38 -b SRR1551047.bam -c target.bed"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :r:b:p:n:f:t:c:uh opt
do
    case $opt in
        b) sbam=$OPTARG;;
        p) pair_id=$OPTARG;;
	r) index_path=$OPTARG;;
	c) capture=$OPTARG;;
	u) umi='umi';;
        h) usage;;
    esac
done

shift $(($OPTIND -1))
baseDir="`dirname \"$0\"`"

if [[ -z $SLURM_CPUS_ON_NODE ]]
then
    SLURM_CPUS_ON_NODE=1
fi

source /etc/profile.d/modules.sh
module load bedtools/2.26.0

if [[ -z $index_path ]] 
then
    index_path='/project/shared/bicf_workflow_ref/human/GRCh38'
fi

if [[ $capture == "${index_path}/clinseq_prj/UTSWV2_2.panelplus.bed" ]]
then
    bedtools coverage -sorted -g  ${index_path}/genomefile.txt -a ${capture} -b ${sbam} -hist > covhist.txt
    grep ^all covhist.txt > genomecov.txt
    sumdepth=`awk '{ sum+= $2*$3;} END {print sum;}' genomecov.txt`
    total=`head -n 1 genomecov.txt |cut -f 4`
    avgdepth=$((${sumdepth}/${total}))
    normals="${index_path}/clinseq_prj/panelofnormals.panel1385V2_2.cnn"
    targets="${index_path}/clinseq_prj/panel1385V2-2.cnvkit_"
    if [[ "$avgdepth" -lt 1000 ]]
    then
	normals="${index_path}/clinseq_prj/panelofnormals.panel1385V2_2.lowcov.cnn"
    fi
elif [[ $capture == "${index_path}/clinseq_prj/UTSWV2.bed" ]]
then 
    normals="${index_path}/clinseq_prj/UTSWV2.normals.cnn"
    targets="${index_path}/clinseq_prj/UTSWV2.cnvkit_"
    if [[ $umi == 'umi' ]]
    then
	normals="${index_path}/clinseq_prj/UTSWV2.uminormals.cnn"
    fi
elif [[ $capture == "${index_path}/clinseq_prj/hemepanelV3.bed" ]]
then
    normals="${index_path}/clinseq_prj/hemepanelV3.panelofnormals.cnn"
    targets="${index_path}/clinseq_prj/hemepanelV3.cnvkit_"
fi	

#echo "$baseDir/../process_scripts/variants/cnvkit.sh -c $capture -b $sbam -p $pair_id -n $normals -t $targets"

bash $baseDir/../process_scripts/variants/cnvkit.sh -c $capture -b $sbam -p $pair_id -n $normals -t $targets
