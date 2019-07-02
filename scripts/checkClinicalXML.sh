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
        c) caseId=$OPTARG;;
	h) usage;;
    esac
done

shift $(($OPTIND -1))
baseDir="`dirname \"$0\"`"

i=$caseId

if [[ ! -f $i.xml ]]
then
    python $baseDir/../IntellispaceDemographics/gatherdemographics.py -i $i -u phg_workflow -p $password -o ${i}.xml
    missing=`grep '><\|D64.9\|N\/A' ${i}.xml`
    if [[ -n $missing ]]
    then
	SUBJECT="SECUREMAIL: Case Missing Data"
	TO="erika.villa@utsouthwestern.edu,Hui.Zheng@UTSouthwestern.edu,Yan.Xu@UTSouthwestern.edu"
	email=$codedir"/bioinformatics_email.txt"
	cat $email | sed 's/Unspecified/'$i'/' | /bin/mail -s "$SUBJECT" "$TO"
    fi
fi
