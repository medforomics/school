#!/bin/bash
#create_directories.sh

usage() {
  echo "-h Help documentation for gatkrunner.sh"
  echo "-p  --projectID"
  echo "-c  --Processing Directory default /project/PHG/PHG_Clinical/processing"
  echo "Example: bash create_directories.sh -p prefix -c /path/processing"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :r:p:c::xh opt
do
    case $opt in
        p) prjid=$OPTARG;;
        c) prodir=$OPTARG;;
	h) usage;;
    esac
done

shift $(($OPTIND -1))
baseDir="`dirname \"$0\"`"

if [[ -z $prodir ]]
then
    prodir="/project/PHG/PHG_Clinical/processing"
fi

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
