#!/bin/bash

set -e
caseID=$1
baseDir="`dirname \"$0\"`"

module load samtools/gcc/1.8

if [[ ! -f ${caseID}.json ]]
then
    curl -o ${caseID}.json "https://nuclia.biohpc.swmed.edu/getSampleAndRunIds?token=$nucliatoken&projectId=$caseID"
fi

jsonout=$(perl ${baseDir}/parse_getsampleJson.pl ${caseID}.json 2>&1)
myarray=($jsonout)
runid=${myarray[2]}
monyear="20${runid:0:4}"

if [[ ! -f  /archive/PHG/PHG_Clinical/toarchive/backups/${monyear}/${myarray[1]}.tar.gz ]]
then
    tar cf /archive/PHG/PHG_Clinical/toarchive/backups/${monyear}/${myarray[1]}.tar /archive/PHG/PHG_Clinical/cases/${myarray[1]}
    gzip /archive/PHG/PHG_Clinical/toarchive/backups/${monyear}/${myarray[1]}.tar 
fi 
if [[ ! -d /archive/PHG/PHG_Clinical/casesAzure/${myarray[1]} ]]
then
    mkdir /archive/PHG/PHG_Clinical/casesAzure/${myarray[1]}
fi
if [[ ! -d /archive/PHG/PHG_Clinical/casesAzure/${myarray[1]}/${myarray[3]} ]]
then
    mkdir /archive/PHG/PHG_Clinical/casesAzure/${myarray[1]}/${myarray[3]}
fi

if [[ -f /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[3]}/${myarray[3]}.consensus.bam ]]
then
    cp /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[3]}/${myarray[3]}.consensus.bam /archive/PHG/PHG_Clinical/casesAzure/${myarray[1]}/${myarray[3]}/${myarray[3]}.bam
    samtools index -@ 4 /archive/PHG/PHG_Clinical/casesAzure/${myarray[1]}/${myarray[3]}/${myarray[3]}.bam
else
    cp /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[3]}/${myarray[3]}.final.bam /archive/PHG/PHG_Clinical/casesAzure/${myarray[1]}/${myarray[3]}/${myarray[3]}.bam
    cp /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[3]}/${myarray[3]}.final.bai /archive/PHG/PHG_Clinical/casesAzure/${myarray[1]}/${myarray[3]}/${myarray[3]}.bam.bai
fi

cp /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[3]}/${myarray[3]}*answer* /archive/PHG/PHG_Clinical/casesAzure/${myarray[1]}/${myarray[3]}
cp /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[1]}.vcf.gz /archive/PHG/PHG_Clinical/casesAzure/${myarray[1]}
cp /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[1]}.TMB.csv /archive/PHG/PHG_Clinical/casesAzure/${myarray[1]}

if [[ -d /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[5]} ]]
then
    cp /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[1]}.concordance.txt /archive/PHG/PHG_Clinical/casesAzure/${myarray[1]}
    cp /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[1]}.utswpass.somatic.vcf* /archive/PHG/PHG_Clinical/casesAzure/${myarray[1]}
    if [[ -f /archive/PHG/PHG_Clinical/casesAzure/${myarray[1]}/${myarray[1]}.utswpass.somatic.vcf.gz ]]
    then
	gunzip /archive/PHG/PHG_Clinical/casesAzure/${myarray[1]}/${myarray[1]}.utswpass.somatic.vcf.gz
    fi
    if [[ ! -d /archive/PHG/PHG_Clinical/casesAzure/${myarray[1]}/${myarray[5]} ]]
    then
	mkdir /archive/PHG/PHG_Clinical/casesAzure/${myarray[1]}/${myarray[5]}
    fi
    if [[ -f /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[5]}/${myarray[5]}.consensus.bam ]]
    then
	samtools index -@ 4 /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[5]}/${myarray[5]}.consensus.bam
	cp /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[5]}/${myarray[5]}.consensus.bam* /archive/PHG/PHG_Clinical/casesAzure/${myarray[1]}/${myarray[5]}
    else
	cp /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[5]}/${myarray[5]}.final.ba* /archive/PHG/PHG_Clinical/casesAzure/${myarray[1]}/${myarray[5]}
    fi	   
fi
    
if [[ -f /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[7]}/${myarray[7]}.bam ]]
then
    if [[ ! -d /archive/PHG/PHG_Clinical/casesAzure/${myarray[1]}/${myarray[7]} ]]
    then
	mkdir /archive/PHG/PHG_Clinical/casesAzure/${myarray[1]}/${myarray[7]}
    fi
    samtools index -@ 4 /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[7]}/${myarray[7]}.bam
    cp /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[7]}/${myarray[7]}.bam* /archive/PHG/PHG_Clinical/casesAzure/${myarray[1]}/${myarray[7]}
    cp /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[7]}/${myarray[7]}.cts /archive/PHG/PHG_Clinical/casesAzure/${myarray[1]}/${myarray[7]}
    cp /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[7]}/${myarray[7]}.fpkm.txt /archive/PHG/PHG_Clinical/casesAzure/${myarray[1]}/${myarray[7]}
    cp /archive/PHG/PHG_Clinical/cases/${myarray[1]}/${myarray[7]}/${myarray[7]}.translocations.answer.txt /archive/PHG/PHG_Clinical/casesAzure/${myarray[1]}/${myarray[7]}
fi

mv /archive/PHG/PHG_Clinical/cases/${myarray[1]} /archive/PHG/PHG_Clinical/toarchive/casesDirs/${myarray[1]}
mv /archive/PHG/PHG_Clinical/casesAzure/${myarray[1]} /archive/PHG/PHG_Clinical/cases
