#!/usr/bin/env nextflow
	
params.input = './fastq'
params.output = '.'
snpeff_vers = 'GRCh38.86';

params.genome="/project/shared/bicf_workflow_ref/human/grch38_cloud/rnaref/"
params.markdups="skip"
params.stranded="0"
params.pairs="pe"
params.align = 'hisat'
params.bamct = "detect"
params.version = 'v5'

params.caseid = 'test'
params.seqrunid = 'runtest'
params.tumorid = 'Tumor'

index_path = file(params.genome)
index_name = "genome"
reffa=file("$params.genome/genome.fa")
gtf_file = file("$params.genome/gencode.gtf")
genenames = file("$params.genome/genenames.txt")

glist=''
if (params.glist) {
   ponopt="-f $params.glist"
}

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

process rtrim {
  executor 'local'
  errorStrategy 'ignore'
  input:
  set pair_id, file(fqs) from reads
  output:
  set env(caseid),pair_id, file("${pair_id}.bam"),file("${pair_id}.bam.bai") into ctbams
  set env(caseid),pair_id, file("${pair_id}.bam") into abundbam
  set env(caseid),pair_id, file("${pair_id}.bam"),file("${pair_id}.bam.bai") into fbbam
  set env(caseid),pair_id, file("${pair_id}.bam"),file("${pair_id}.alignerout.txt") into qcbam
  set env(caseid), pair_id, file("${pair_id}.trim.R1.fastq.gz"),file("${pair_id}.trim.R2.fastq.gz") into fusionfq
  script:
  """
  caseid=${params.caseid}
  bash $baseDir/process_scripts/preproc_fastq/trimgalore.sh -p ${pair_id} ${fqs}
  bash $baseDir/process_scripts/alignment/rnaseqalign.sh -a $params.align -p $pair_id -r $index_path -x ${pair_id}.trim.R1.fastq.gz -y ${pair_id}.trim.R2.fastq.gz
  """
}
process starfusion {
  errorStrategy 'ignore'
  publishDir "$params.output/$caseid/$pair_id", mode: 'copy'
  input:
  set caseid,pair_id, file(fq1), file(fq2) from fusionfq
  output:
  file("${pair_id}*txt") into fusionout
  script:
  """
  bash $baseDir/process_scripts/alignment/starfusion.sh -p ${pair_id} -r ${index_path} -a ${fq1} -b ${fq2} -m trinity -f
  """
}
process bamct {
  executor 'local'
  errorStrategy 'ignore'
  publishDir "$params.output/$caseid/$pair_id", mode: 'copy'
  input:
  set caseid,pair_id, file(rbam),file(ridx) from ctbams
  output:
  file("${pair_id}.bamreadct.txt") into ctreads
  when:
  params.bamct == "detect"
  script:
  """
  /project/shared/bicf_workflow_ref/seqprg/bam-readcount/bin/bam-readcount -w 0 -q 0 -b 25 -f ${index_path}/genome.fa $rbam > ${pair_id}.bamreadct.txt
  """
}
process alignqc {
  executor 'local'
  errorStrategy 'ignore'
  publishDir "$params.output/$caseid/$pair_id", mode: 'copy'

  input:
  set caseid,pair_id, file(bam), file(hsout) from qcbam
  
  output:
  set file("${pair_id}_fastqc.zip"),file("${pair_id}_fastqc.html") into fastqc
  file("${pair_id}.sequence.stats.txt") into alignstats
  script:
  """
  bash $baseDir/process_scripts/alignment/bamqc.sh -p ${pair_id} -b ${bam} -n rna
  perl $baseDir/scripts/sequenceqc_rna.pl -r ${index_path} -e ${params.version} *.flagstat.txt
  """
}
process geneabund {
  errorStrategy 'ignore'
  executor 'local'
  publishDir "$params.output/$caseid/$pair_id", mode: 'copy'
  input:
  set caseid,pair_id, file(sbam) from abundbam
  output:
  file("${pair_id}.cts")  into counts
  file("${pair_id}_stringtie") into strcts
  file("${pair_id}.fpkm.txt") into fpkm
  """
  source /etc/profile.d/modules.sh
  bash $baseDir/process_scripts/diff_exp/geneabundance.sh -s $params.stranded -g ${gtf_file} -p ${pair_id} -b ${sbam} $glist -f 1
  """
}
process fb {
  queue '32GB'
  errorStrategy 'ignore'
  publishDir "$params.output/$caseid", mode: 'copy'

  input:
  set caseid,$pair_id,file(gbam),file(gidx) from fbbam
  output:
  set caseid,file("${caseid}_${params.seqrunid}.rna.vcf.gz") into fbvcf
  script:
  """
  bash $baseDir/process_scripts/variants/germline_vc.sh -r $index_path -p $caseid -a fb
  bash $baseDir/process_scripts/variants/uni_norm_annot.sh -g $snpeff_vers -r $index_path -p ${caseid}.fb -v ${caseid}.fb.vcf.gz 
  mv ${caseid}.fb.vcf.gz ${caseid}_${params.seqrunid}.rna.vcf.gz
  """
}
