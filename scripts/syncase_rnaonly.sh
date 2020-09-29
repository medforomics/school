#!/bin/bash

set -e
caseID=$1
nucliatoken=$2
baseDir="`dirname \"$0\"`"

source ${baseDir}/../azure_credentials
module load samtools/gcc/1.8

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

mkdir -p /archive/PHG/PHG_Clinical/casesTemp/${myarray[1]}

if [[ -f /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[1]}.vcf.gz ]]
then
    cp /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[1]}.vcf.gz /archive/PHG/PHG_Clinical/casesTemp/${myarray[1]}
   cp /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[1]}.maf /archive/PHG/PHG_Clinical/casesTemp/${myarray[1]}
   cp /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[1]}.TMB.csv /archive/PHG/PHG_Clinical/casesTemp/${myarray[1]}
fi
   
if [[ -f /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[1]}.mutational_signature.txt ]]
then
    cp /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[1]}.mutational_signature.txt /archive/PHG/PHG_Clinical/casesTemp/${myarray[1]}
    cp /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[1]}.mutational_signature.png /archive/PHG/PHG_Clinical/casesTemp/${myarray[1]}
fi 

if [[ -d /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[3]} ]]
then
    mkdir -p /archive/PHG/PHG_Clinical/casesTemp/${myarray[1]}/${myarray[3]}
    cp /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[3]}/${myarray[3]}*answer* /archive/PHG/PHG_Clinical/casesTemp/${myarray[1]}/${myarray[3]}
fi

if [[ -f /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[3]}/${myarray[3]}.consensus.bam ]]
then
    cp /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[3]}/${myarray[3]}.consensus.bam /archive/PHG/PHG_Clinical/casesTemp/${myarray[1]}/${myarray[3]}/${myarray[3]}.bam
    samtools index -@ 4 /archive/PHG/PHG_Clinical/casesTemp/${myarray[1]}/${myarray[3]}/${myarray[3]}.bam
elif [[ -f /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[3]}/${myarray[3]}.final.bam ]]
then
    cp /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[3]}/${myarray[3]}.final.bam /archive/PHG/PHG_Clinical/casesTemp/${myarray[1]}/${myarray[3]}/${myarray[3]}.bam
    cp /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[3]}/${myarray[3]}.final.bai /archive/PHG/PHG_Clinical/casesTemp/${myarray[1]}/${myarray[3]}/${myarray[3]}.bam.bai
fi
if [[ -d /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[5]} ]]
then
    cp /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[1]}.concordance.txt /archive/PHG/PHG_Clinical/casesTemp/${myarray[1]}
    cp /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[1]}.utswpass.somatic.vcf* /archive/PHG/PHG_Clinical/casesTemp/${myarray[1]}
    if [[ -f /archive/PHG/PHG_Clinical/casesTemp/${myarray[1]}/${myarray[1]}.utswpass.somatic.vcf.gz ]]
    then
	gunzip /archive/PHG/PHG_Clinical/casesTemp/${myarray[1]}/${myarray[1]}.utswpass.somatic.vcf.gz
    fi
    mkdir -p /archive/PHG/PHG_Clinical/casesTemp/${myarray[1]}/${myarray[5]}
    if [[ -f /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[5]}/${myarray[5]}.consensus.bam ]]
    then
	cp /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[5]}/${myarray[5]}.consensus.bam /archive/PHG/PHG_Clinical/casesTemp/${myarray[1]}/${myarray[5]}/${myarray[5]}.bam
	samtools index -@ 4 /archive/PHG/PHG_Clinical/casesTemp/${myarray[1]}/${myarray[5]}/${myarray[5]}.bam
    else
	cp /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[5]}/${myarray[5]}.final.bam /archive/PHG/PHG_Clinical/casesTemp/${myarray[1]}/${myarray[5]}/${myarray[5]}.bam
	samtools index -@ 4 /archive/PHG/PHG_Clinical/casesTemp/${myarray[1]}/${myarray[5]}/${myarray[5]}.bam
    fi	   
fi

if [[ -f /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[7]}/${myarray[7]}.bam ]]
then
    mkdir -p /archive/PHG/PHG_Clinical/casesTemp/${myarray[1]}/${myarray[7]}
    if [[ -f /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[7]}/${myarray[7]}.bamreadct.txt ]]
    then
	gzip /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[7]}/${myarray[7]}.bamreadct.txt
    fi
    cp /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[7]}/${myarray[7]}.bam /archive/PHG/PHG_Clinical/casesTemp/${myarray[1]}/${myarray[7]}
    samtools index -@ 4 /archive/PHG/PHG_Clinical/casesTemp/${myarray[1]}/${myarray[7]}/${myarray[7]}.bam
    cp /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[7]}/${myarray[7]}.cts /archive/PHG/PHG_Clinical/casesTemp/${myarray[1]}/${myarray[7]}
    cp /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[7]}/${myarray[7]}.fpkm.txt /archive/PHG/PHG_Clinical/casesTemp/${myarray[1]}/${myarray[7]}
    cp /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[7]}/${myarray[7]}.translocations.answer.txt /archive/PHG/PHG_Clinical/casesTemp/${myarray[1]}/${myarray[7]}
fi

mv /archive/PHG/PHG_Clinical/cases/${myarray[1]} /archive/PHG/PHG_Clinical/toarchive/caseDirs/${myarray[1]}
mv /archive/PHG/PHG_Clinical/casesTemp/${myarray[1]} /archive/PHG/PHG_Clinical/cases
rsync -avz /archive/PHG/PHG_Clinical/cases/${myarray[1]} answerbe@198.215.54.71:/swnas/cases

cd /archive/PHG/PHG_Clinical/cases
az storage blob upload-batch -d cases -s /archive/PHG/PHG_Clinical/cases/${myarray[1]} --destination-path ${year}/${myarray[1]} > /archive/PHG/PHG_Clinical/toarchive/transfer_logs/${year}/${seqrunid}/${caseID}.log
