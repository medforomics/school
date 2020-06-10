#!/usr/bin/env nextflow

params.input = './fastq'
params.output = './analysis'

def somatic = [:]
def ids = []
def reads = []

if ( params.caseid ) {
   params.tumorid = 'Tumor'
   somatic[params.caseid] = false		
   params.tfq = "${params.input}/${params.tumorid}.R{1,2}.fastq.gz"
   reads << tuple(params.caseid,params.tumorid,file(params.tfq))
   if( params.normalid ) {
      somatic[params.caseid] = true
      params.nfq = "${params.input}/${params.normalid}.R{1,2}.fastq.gz"
      ids << tuple(params.caseid,params.tumorid,params.normalid)
      ids << tuple(params.caseid,params.tumorid,params.normalid)
      reads << tuple(params.caseid,params.normalid,file(params.nfq))
   } else {
      ids << tuple(params.caseid,params.tumorid,'')
   }
} else {
   params.design="$params.input/design.txt"
   params.fastqs="$params.input/*.fastq.gz"
   fastqs=file(params.fastqs)
   design_file=file(params.design)
   def fileMap = [:]
   fastqs.each {
      final fileName = it.getFileName().toString()
      prefix = fileName.lastIndexOf('/')
      fileMap[fileName] = it
   }
   new File(params.design).withReader { reader ->
      def hline = reader.readLine()
      def header = hline.split("\t")
      cidx = header.findIndexOf{it == 'CaseID'};
      tidx = header.findIndexOf{it == 'TumorID'};
      nidx = header.findIndexOf{it == 'NormalID'};
      fidx = header.findIndexOf{it == 'SampleID'};
      oneidx = header.findIndexOf{it == 'FQR1'};
      twoidx = header.findIndexOf{it == 'FQR2'};
      while (line = reader.readLine()) {
    	   def row = line.split("\t")
      if (fileMap.get(row[oneidx]) != null) {
       	   somatic[row[cidx]] = true
      	   if (row[nidx] == '') {
       	       somatic[row[cidx]] = false
       	   }
     	   reads << tuple(row[cidx],row[fidx],[fileMap.get(row[oneidx]),fileMap.get(row[twoidx])])
     	   ids << tuple(row[cidx],row[tidx],row[nidx])
      }	
   }
 }
}
if( ! reads) { error "Didn't match any input files with entries in the design file" }

reads.println()
ids.println()

process dtrim_align {
   executor 'local'
   echo true
   errorStrategy 'ignore'
   publishDir "$params.output/$caseid/$pair_id", mode: 'copy'
   input:
   set caseid,pair_id,file(fqs) from reads
   set caseid,tumorid,normalid from ids
   output:     
   set caseid,tumorid,normalid,file("${pair_id}.txt"), file("${pair_id}.twice.txt") into oribam
   script:  
   if ( somatic[caseid] == true )
   """
   echo somatic $caseid $normalid $tumorid $pair_id $fqs
   touch ${pair_id}.txt
   zcat $fqs >${pair_id}.twice.txt
   """
   else
   """
   echo tumoronly $caseid 'nonormal' $tumorid $pair_id $fqs
   touch ${pair_id}.txt
   zcat $fqs >${pair_id}.twice.txt
   """
}

oribam
   .groupTuple(by:[0,1,2])
   .set { mutectbam }

process runVC {
  executor 'local'
  echo true
  errorStrategy 'ignore'
  publishDir "$params.output/$caseid/$pair_id", mode: 'copy'
  input:
   set caseid,tid,nid,file(ssbam),file(ssidx) from mutectbam
   output:     
   stdout into result3
   script:  
   """
   echo $caseid $tid $nid $ssbam $ssidx
   """
}
