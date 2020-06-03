nextflow.preview.dsl=2



process trim {
  executor 'local'
  errorStrategy 'ignore'
  input:
  set subjid, pair_id, file(read1), file(read2)
  output:
  set subjid, pair_id, file("${pair_id}.trim.R1.fastq.gz"),file("${pair_id}.trim.R2.fastq.gz"),file("${pair_id}.trimreport.txt")
  script:
  """
  bash $baseDir/process_scripts/preproc_fastq/trimgalore.sh -p ${pair_id} -a ${read1} -b ${read2}
  """
}
