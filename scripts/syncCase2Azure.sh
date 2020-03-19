#!/bin/bash

set -e
caseID=$1
nucliatoken=$2
wrkdir=$3

if [[ -z $wrkdir ]]
then
    wkdir= /archive/PHG/PHG_Clinical
fi

baseDir="`dirname \"$0\"`"

module load samtools/gcc/1.8 azure/2.0.72 
export HTTP_PROXY="http://proxy.swmed.edu:3128"
export HTTPS_PROXY="http://proxy.swmed.edu:3128"
export http_proxy="http://proxy.swmed.edu:3128"
export https_proxy="http://proxy.swmed.edu:3128"
export AZURE_STORAGE_ACCOUNT=swazrstrseq
export AZURE_STORAGE_KEY=SxfAj0wkNDyQVKcP5ChoTDEd7J8bZm2zAqh0vNE2YRAQxFGrTLc1wvlWv0IYrS1p6thpBCvBLGHbVQmiu1/XqQ==

if [[ ! -f ${caseID}.json ]]
then
    curl -o ${caseID}.json "https://nuclia.biohpc.swmed.edu/getSampleAndRunIds?token=$nucliatoken&projectId=$caseID"
fi

jsonout=$(perl ${baseDir}/parse_getsampleJson.pl ${caseID}.json 2>&1)
myarray=($jsonout)
if [[ ${myarray[2]} != "NA" ]]
then
    runid=${myarray[2]}
elif [[ ${myarray[6]} != "NA" ]]
then
    runid=${myarray[6]}
fi
monyear="20${runid:0:4}"

if [[ ! -d $wkdir/casesTemp/${myarray[1]} ]]
then
    mkdir $wkdir/casesTemp/${myarray[1]}
fi
if [[ -d $wkdir/cases/${myarray[1]}/${myarray[3]} ]]
then
    if [[ ! -d $wkdir/casesTemp/${myarray[1]}/${myarray[3]} ]]
    then
	mkdir $wkdir/casesTemp/${myarray[1]}/${myarray[3]}
    fi
    cp $wkdir/cases/${myarray[1]}/${myarray[3]}/${myarray[3]}*answer* $wkdir/casesTemp/${myarray[1]}/${myarray[3]}
    cp $wkdir/cases/${myarray[1]}/${myarray[3]}/${myarray[3]}.ballelefreq.txt $wkdir/casesTemp/${myarray[1]}/${myarray[3]}
    cp $wkdir/cases/${myarray[1]}/${myarray[1]}.vcf.gz $wkdir/casesTemp/${myarray[1]}
    if -f [[ $wkdir/cases/${myarray[1]}/${myarray[1]}.translocations.answer.txt ]]
    then
	cp $wkdir/cases/${myarray[1]}/${myarray[1]}.translocations.answer.txt $wkdir/casesTemp/${myarray[1]}
    fi
fi
if [[ -f $wkdir/cases/${myarray[1]}/${myarray[3]}/${myarray[3]}.consensus.bam ]]
then
    cp $wkdir/cases/${myarray[1]}/${myarray[3]}/${myarray[3]}.consensus.bam $wkdir/casesTemp/${myarray[1]}/${myarray[3]}/${myarray[3]}.bam
    samtools index -@ 4 $wkdir/casesTemp/${myarray[1]}/${myarray[3]}/${myarray[3]}.bam
elif [[ -f $wkdir/cases/${myarray[1]}/${myarray[3]}/${myarray[3]}.final.bam ]]
then
    cp $wkdir/cases/${myarray[1]}/${myarray[3]}/${myarray[3]}.final.bam $wkdir/casesTemp/${myarray[1]}/${myarray[3]}/${myarray[3]}.bam
    cp $wkdir/cases/${myarray[1]}/${myarray[3]}/${myarray[3]}.final.bai $wkdir/casesTemp/${myarray[1]}/${myarray[3]}/${myarray[3]}.bam.bai
fi
if [[ -d $wkdir/cases/${myarray[1]}/${myarray[5]} ]]
then
    cp $wkdir/cases/${myarray[1]}/${myarray[1]}.concordance.txt $wkdir/casesTemp/${myarray[1]}
    cp $wkdir/cases/${myarray[1]}/${myarray[1]}.utswpass.somatic.vcf* $wkdir/casesTemp/${myarray[1]}
    cp $wkdir/cases/${myarray[1]}/${myarray[1]}.TMB.csv $wkdir/casesTemp/${myarray[1]}
    if [[ -f $wkdir/casesTemp/${myarray[1]}/${myarray[1]}.utswpass.somatic.vcf.gz ]]
    then
	gunzip $wkdir/casesTemp/${myarray[1]}/${myarray[1]}.utswpass.somatic.vcf.gz
    fi
    if [[ ! -d $wkdir/casesTemp/${myarray[1]}/${myarray[5]} ]]
    then
	mkdir $wkdir/casesTemp/${myarray[1]}/${myarray[5]}
    fi
    if [[ -f $wkdir/cases/${myarray[1]}/${myarray[5]}/${myarray[5]}.consensus.bam ]]
    then
	cp $wkdir/cases/${myarray[1]}/${myarray[5]}/${myarray[5]}.consensus.bam $wkdir/casesTemp/${myarray[1]}/${myarray[5]}/${myarray[5]}.bam
	samtools index -@ 4 $wkdir/casesTemp/${myarray[1]}/${myarray[5]}/${myarray[5]}.bam
    else
	cp $wkdir/cases/${myarray[1]}/${myarray[5]}/${myarray[5]}.final.bam $wkdir/casesTemp/${myarray[1]}/${myarray[5]}/${myarray[5]}.bam
	samtools index -@ 4 $wkdir/casesTemp/${myarray[1]}/${myarray[5]}/${myarray[5]}.bam
    fi	   
fi
    
if [[ -f $wkdir/cases/${myarray[1]}/${myarray[7]}/${myarray[7]}.bam ]]
then
    if [[ ! -d $wkdir/casesTemp/${myarray[1]}/${myarray[7]} ]]
    then
	mkdir $wkdir/casesTemp/${myarray[1]}/${myarray[7]}
    fi
    if [[ -f $wkdir/cases/${myarray[1]}/${myarray[7]}/${myarray[7]}.bamreadct.txt ]]
    then
	gzip $wkdir/cases/${myarray[1]}/${myarray[7]}/${myarray[7]}.bamreadct.txt
    fi
    cp $wkdir/cases/${myarray[1]}/${myarray[7]}/${myarray[7]}.bam $wkdir/casesTemp/${myarray[1]}/${myarray[7]}
    samtools index -@ 4 $wkdir/casesTemp/${myarray[1]}/${myarray[7]}/${myarray[7]}.bam
    cp $wkdir/cases/${myarray[1]}/${myarray[7]}/${myarray[7]}.cts $wkdir/casesTemp/${myarray[1]}/${myarray[7]}
    cp $wkdir/cases/${myarray[1]}/${myarray[7]}/${myarray[7]}.fpkm.txt $wkdir/casesTemp/${myarray[1]}/${myarray[7]}
    if [[ -f $wkdir/cases/${myarray[1]}/${myarray[1]}.translocations.answer.txt ]]
    then
	cp $wkdir/cases/${myarray[1]}/${myarray[1]}.translocations.answer.txt $wkdir/casesTemp/${myarray[1]}/${myarray[7]}/${myarray[7]}.translocations.answer.txt
    else
        cp $wkdir/cases/${myarray[1]}/${myarray[7]}/${myarray[7]}.translocations.answer.txt $wkdir/casesTemp/${myarray[1]}/${myarray[7]}
    fi
    
fi

mv $wkdir/cases/${myarray[1]} $wkdir/toarchive/caseDirs/${myarray[1]}
mv $wkdir/casesTemp/${myarray[1]} $wkdir/cases
rsync -avz $wkdir/cases/${myarray[1]} answerbe@198.215.54.71:/swnas/cases

if [[ ! -f  $wkdir/toarchive/backups/${monyear}/${myarray[1]}.tar.gz ]]
then
    tar cf $wkdir/toarchive/backups/${monyear}/${myarray[1]}.tar $wkdir/toarchive/caseDirs/${myarray[1]}
    gzip $wkdir/toarchive/backups/${monyear}/${myarray[1]}.tar 
fi 

isfalse=`az storage container exists -n archive${monyear} | grep false`
if [[ -n $isfalse ]]
then
    az storage container create -n archive${monyear} --fail-on-exist
fi

az storage blob upload -c archive${monyear} -f $wkdir/toarchive/backups/${monyear}/${myarray[1]}.tar.gz -n ${myarray[1]}.tar.gz
