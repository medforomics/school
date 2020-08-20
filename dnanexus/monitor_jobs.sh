#!/bin/bash
#monitor_jobs.sh

set -e
baseDir="`dirname \"$0\"`"
module load dxtoolkit/python27/0.294.0 
jobbase='/project/PHG/PHG_Clinical/cloud'
processing='/project/PHG/PHG_Clinical/devel/processing'

for i in ${jobbase}/pending/*.joblist.txt
do
    filename="$(basename -- $i)"
    runid=$(echo $filename | cut -f 1 -d  '.')
    caseid=$(echo $filename | cut -f 2 -d  '.')
    complete=0
    fail=0
    numjobs=$(grep -c analysis $i)
    mkdir -p ${processing}/${runid}/${caseid}
    cd ${processing}/${runid}/${caseid}
    while read j
    do
	dx describe $j --json > $j.json
	state=$(cat $j.json | jq -r '.state')
	if [[ $state == 'done' ]]
	then
	    complete=$((complete+1))
	elif [[ $state == 'failed' ]]
	then
	    fail=1
	    SUBJECT="SECUREMAIL: Case Job Failed: $j"
            TO="ngsclialab@UTSouthwestern.edu,Chelsea.Raulerson@UTSouthwestern.edu"
            email="${baseDir}/dnanexus_email.txt"
	    #cat $email | sed 's/Unspecified/'$j'/' 
	    cat $email | sed 's/Unspecified/'$j'/' | /bin/mail -s "$SUBJECT" "$TO"
	fi
    done<$i
    if [[ $complete == $numjobs ]]
    then
	mv $i $jobbase/complete
	cd ${processing}/${runid}/${caseid}
	dx ls /${runid}/${caseid} --folders > folder.list.txt
	while read outf
	do
	    dx download /${runid}/${caseid}/$outf -r
	    cd $outf
	    ls *.tar.gz |awk '{print "tar xvfz",$1}' |sh
	    cd ..
	done<folder.list.txt 
    elif [[ $fail == 1 ]]
    then
	mv $i $jobbase/archive
    fi
done

#### Code to get files and exclude certain files
#dx find data --path /${runid}/${caseid} --json > files.json
#filemap=$(cat files.json | jq -c '.[] | [.describe.name, .describe.id]')
#for i in $filemap
#do
#line=$(echo $i | column -t -s'[],"')
#myarray=($line)
#if [[ ! ${myarray[0]} =~ 'fastq.gz' ]] && [[ ! ${myarray[0]} =~ 'final.ba' ]]
#then
#echo "dx download ${myarray[1]}"
#else
#echo "skipping ${myarray[0]}"
#fi
#done
