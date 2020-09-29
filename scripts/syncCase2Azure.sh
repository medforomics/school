#!/bin/bash

set -e
caseID=$1
nucliatoken=$2
wkdir=$3
testing=$4

if [[ -z $wkdir ]]
then
    wkdir=/archive/PHG/PHG_Clinical
fi
baseDir="`dirname \"$0\"`"

module load samtools/gcc/1.8

if [[ ! -f ${caseID}.json ]]
then
    curl -o ${caseID}.json "https://nuclia.biohpc.swmed.edu/getSampleAndRunIds?token=$nucliatoken&projectId=$caseID"
fi

jsonout=$(perl ${baseDir}/parse_getsampleJson.pl ${caseID}.json 2>&1)
myarray=($jsonout)
subject=${myarray[1]}
tumor_id=${myarray[3]}
dna_runid=${myarray[2]}
rna_runid=${myarray[6]}
normal_id=${myarray[5]}
rnaseq_id=${myarray[7]}

if [[ ${myarray[2]} != "NA" ]]
then
    runid=${myarray[2]}
elif [[ ${myarray[6]} != "NA" ]]
then
    runid=${myarray[6]}
fi
monyear="20${runid:0:4}"

mkdir -p ${wkdir}/casesTemp/${subject}
cp ${wkdir}/cases/${subject}/${subject}.vcf.gz ${wkdir}/casesTemp/${subject}
cp ${wkdir}/cases/${subject}/${subject}.maf ${wkdir}/casesTemp/${subject}
cp ${wkdir}/cases/${subject}/${subject}.viral_results.txt ${wkdir}/casesTemp/${subject}
cp ${wkdir}/cases/${subject}/${subject}.TMB.csv ${wkdir}/casesTemp/${subject}

if [[ -f /archive/PHG/PHG_Clinical/cases/${subject}/${subject}.mutational_signature.txt ]]
then
    cp /archive/PHG/PHG_Clinical/cases/${subject}/${subject}.mutational_signature.txt /archive/PHG/PHG_Clinical/casesTemp/${subject}
    cp /archive/PHG/PHG_Clinical/cases/${subject}/${subject}.mutational_signature.png /archive/PHG/PHG_Clinical/casesTemp/${subject}
fi 

if [[ -d ${wkdir}/cases/${subject}/${tumor_id} ]]
then
    mkdir -p ${wkdir}/casesTemp/${subject}/${tumor_id}
    cp ${wkdir}/cases/${subject}/${tumor_id}/${tumor_id}*answer* ${wkdir}/casesTemp/${subject}/${tumor_id}
    cp ${wkdir}/cases/${subject}/${tumor_id}/${tumor_id}.ballelefreq.txt ${wkdir}/casesTemp/${subject}/${tumor_id}
    if [[ -f ${wkdir}/cases/${subject}/${subject}.translocations.answer.txt ]]
    then
	cp ${wkdir}/cases/${subject}/${subject}.translocations.answer.txt ${wkdir}/casesTemp/${subject}
    fi
fi
if [[ -f ${wkdir}/cases/${subject}/${tumor_id}/${tumor_id}.consensus.bam ]]
then
    cp ${wkdir}/cases/${subject}/${tumor_id}/${tumor_id}.consensus.bam ${wkdir}/casesTemp/${subject}/${tumor_id}/${tumor_id}.bam
    samtools index -@ 4 ${wkdir}/casesTemp/${subject}/${tumor_id}/${tumor_id}.bam
elif [[ -f ${wkdir}/cases/${subject}/${tumor_id}/${tumor_id}.final.bam ]]
then
    cp ${wkdir}/cases/${subject}/${tumor_id}/${tumor_id}.final.bam ${wkdir}/casesTemp/${subject}/${tumor_id}/${tumor_id}.bam
    cp ${wkdir}/cases/${subject}/${tumor_id}/${tumor_id}.final.bai ${wkdir}/casesTemp/${subject}/${tumor_id}/${tumor_id}.bam.bai
fi
if [[ -d ${wkdir}/cases/${subject}/${normal_id} ]]
then
    cp ${wkdir}/cases/${subject}/${subject}.concordance.txt ${wkdir}/casesTemp/${subject}
    cp ${wkdir}/cases/${subject}/${subject}.utswpass.somatic.vcf* ${wkdir}/casesTemp/${subject}
    if [[ -f ${wkdir}/casesTemp/${subject}/${subject}.utswpass.somatic.vcf.gz ]]
    then
	gunzip ${wkdir}/casesTemp/${subject}/${subject}.utswpass.somatic.vcf.gz
    fi
    mkdir -p ${wkdir}/casesTemp/${subject}/${normal_id}
    if [[ -f ${wkdir}/cases/${subject}/${normal_id}/${normal_id}.consensus.bam ]]
    then
	cp ${wkdir}/cases/${subject}/${normal_id}/${normal_id}.consensus.bam ${wkdir}/casesTemp/${subject}/${normal_id}/${normal_id}.bam
	samtools index -@ 4 ${wkdir}/casesTemp/${subject}/${normal_id}/${normal_id}.bam
    else
	cp ${wkdir}/cases/${subject}/${normal_id}/${normal_id}.final.bam ${wkdir}/casesTemp/${subject}/${normal_id}/${normal_id}.bam
	samtools index -@ 4 ${wkdir}/casesTemp/${subject}/${normal_id}/${normal_id}.bam
    fi	   
fi
if [[ -f ${wkdir}/cases/${subject}/${rnaseq_id}/${rnaseq_id}.bam ]]
then
    mkdir -p ${wkdir}/casesTemp/${subject}/${rnaseq_id}
    if [[ -f ${wkdir}/cases/${subject}/${rnaseq_id}/${rnaseq_id}.bamreadct.txt ]]
    then
	gzip ${wkdir}/cases/${subject}/${rnaseq_id}/${rnaseq_id}.bamreadct.txt
    fi
    cp ${wkdir}/cases/${subject}/${rnaseq_id}/${rnaseq_id}.bam ${wkdir}/casesTemp/${subject}/${rnaseq_id}
    samtools index -@ 4 ${wkdir}/casesTemp/${subject}/${rnaseq_id}/${rnaseq_id}.bam
    cp ${wkdir}/cases/${subject}/${rnaseq_id}/${rnaseq_id}.cts ${wkdir}/casesTemp/${subject}/${rnaseq_id}
    cp ${wkdir}/cases/${subject}/${rnaseq_id}/${rnaseq_id}.fpkm.txt ${wkdir}/casesTemp/${subject}/${rnaseq_id}
    if [[ -f ${wkdir}/cases/${subject}/${subject}.translocations.answer.txt ]]
    then
	cp ${wkdir}/cases/${subject}/${subject}.translocations.answer.txt ${wkdir}/casesTemp/${subject}/${rnaseq_id}/${rnaseq_id}.translocations.answer.txt
    else
        cp ${wkdir}/cases/${subject}/${rnaseq_id}/${rnaseq_id}.translocations.answer.txt ${wkdir}/casesTemp/${subject}/${rnaseq_id}
    fi
fi

if [[ -z $testing ]]
then
    mv $wkdir/cases/${subject} $wkdir/toarchive/caseDirs/${subject}
    mv $wkdir/casesTemp/${subject} $wkdir/cases
fi

if [[ -z $testing ]]
then
    rsync -avz $wkdir/cases/${subject} answerbe@198.215.54.71:/swnas/cases
    source ${baseDir}/../azure_credentials
    az storage blob upload-batch -d seqruns -s $wkdir/cases/${subject} --destination-path ${year}/${subject} >  /archive/PHG/PHG_Clinical/toarchive/transfer_logs/${year}/${seqrunid}/${caseID}.log
else
    rsync -avz $wkdir/casesTemp/${subject} answerbe@198.215.54.71:/swnas/cases/answer_devel
fi
