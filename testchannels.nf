#!/usr/bin/env nextflow


params.input = './fastq'
params.output = './analysis'

params.caseid = 'test'
params.seqrunid = 'runtest'
params.vcfid = "${params.caseid}.${params.seqrunid}"
params.tumorid = 'Tumor'

params.tfq = "${params.input}/${params.tumorid}.R{1,2}.fastq.gz"
Channel
  .fromPath(params.tfq)
  .ifEmpty { error "Cannot find any reads matching: ${params.tfq}" }
  .set { fqfiles }
  
pnames = Channel
  .from([params.tumorid,params.tumorid])
  .merge( fqfiles )
  .groupTuple(by:0)
  .set { reads }

if( params.normalid ) {
    params.nfq = "${params.input}/${params.normalid}.R{1,2}.fastq.gz"
    Channel
	.fromPath([params.tfq, params.nfq])
	.ifEmpty { error "Cannot find any reads matching: ${params.tfq}" }
    	.set { fqfiles }
    pnames = Channel
      .from([params.tumorid, params.tumorid, params.normalid,params.normalid])
      .merge( fqfiles )
      .groupTuple(by:0)
      .set { reads }

}

process dtrim_align {
  executor 'local'
  errorStrategy 'ignore'
  publishDir "$params.output/$params.caseid/$pair_id", mode: 'copy'
  input:
  set pair_id,file(fqs) from reads
  output: 
  stdout into result3
  
  """
  echo $pair_id ${fqs}
  """
}
result3.println()
