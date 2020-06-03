#!/bin/bash


for i in /project/PHG/PHG_Clinical/illumina/1905*; do
    if [[ -d "$i" ]]
    then
	runid=`echo "$i" | cut -f 6 -d '/'`
	rundate=`echo "$runid" | cut -f 1 -d '_'`
	destdir=${rundate:0:4}
	if [[ -f /archive/PHG/PHG_Clinical/toarchive/backups/20${destdir}/${runid}.tar.gz ]]
	then
	    echo "File exists:  /archive/PHG/PHG_Clinical/toarchive/backups/20${destdir}/${runid}.tar.gz"
	else
	    tar cf /archive/PHG/PHG_Clinical/toarchive/backups/20${destdir}/${runid}.tar $i
	    gzip /archive/PHG/PHG_Clinical/toarchive/backups/20${destdir}/${runid}.tar
	fi
    fi
done
