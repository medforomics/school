#!/usr/bin/env nextflow
	
params.input = './fastq'
params.output = './analysis'
snpeff_vers = 'GRCh38.86';

params.genome="/project/shared/bicf_workflow_ref/human/grch38_cloud/rnaref/"
params.geneinfo="/project/shared/bicf_workflow_ref/human/gene_info.human.txt"

params.stranded="0"
params.pairs="pe"
params.align = 'hisat'
params.bamct = "detect"
params.version = 'v5'

params.seqrunid = 'runtest'

ginfo=file(params.geneinfo)
index_path = file(params.genome)
index_name = "genome"
reffa=file("$params.genome/genome.fa")
gtf_file = file("$params.genome/gencode.gtf")
genenames = file("$params.genome/genenames.txt")

glist=''
if (params.glist) {
   glist="-f $params.glist"
}
umiopt=''
if (params.umi) {
   umiopt=" -u"
}
repoDir=workflow.projectDir
if (params.repoDir) {
   repoDir=params.repoDir
}

def reads = []

if (params.caseid) {
   params.tumorid = 'Tumor'
   params.tfq = "${params.input}/${params.tumorid}.R{1,2}.fastq.gz"
   reads << tuple(params.caseid,params.tumorid,file(params.tfq))
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
       fidx = header.findIndexOf{it == 'SampleID'};
       oneidx = header.findIndexOf{it == 'FqR1'};
       twoidx = header.findIndexOf{it == 'FqR2'};
       while (line = reader.readLine()) {
  	   def row = line.split("\t")
      	   if (fileMap.get(row[oneidx]) != null) {
     	      reads << tuple(row[cidx],row[fidx],[fileMap.get(row[oneidx]),fileMap.get(row[twoidx])])
      	   }	
       }
   }
}

if( ! reads) { error "Didn't match any input files with entries in the design file" }

process rtrim {
  errorStrategy 'ignore'
  label 'trim'
  publishDir "$params.output/$caseid/rnaout", mode: 'copy'

  input:
  set caseid,sampleid, file(fqs) from reads
  output:
  set caseid,sampleid,file("${sampleid}.trim.R1.fastq.gz"),file("${sampleid}.trim.R2.fastq.gz") into fusionfq
  set caseid,sampleid,file("${sampleid}.trim.R*.fastq.gz") into treads

  script:
  """
  bash ${repoDir}/process_scripts/preproc_fastq/trimgalore.sh -p ${sampleid} ${fqs}
  """
}
process ralign {
  errorStrategy 'ignore'
  label 'ralign'
  publishDir "$params.output/$caseid/rnaout", mode: 'copy'
  input:
  set caseid,sampleid, file(fqs) from treads
  output:
  set caseid,sampleid,file("${sampleid}.bam") into abundbam
  set caseid,sampleid,file("${sampleid}.bam"),file("${sampleid}.bam.bai") into fbbam
  set caseid,sampleid,file("${sampleid}.bam"),file("${sampleid}.bam.bai") into ctbam
  set caseid,sampleid,file("${sampleid}.bam"),file("${sampleid}.alignerout.txt") into qcbam
  script:
  """
  bash ${repoDir}/process_scripts/alignment/rnaseqalign.sh -a $params.align -p $sampleid -r $index_path $umiopt ${fqs}
  """
}
process starfusion {
  errorStrategy 'ignore'
  label 'starfusion'
  publishDir "$params.output/$caseid/rnaout", mode: 'copy'
  input:
  set caseid,sampleid,file(fq1), file(fq2) from fusionfq
  output:
  file("${sampleid}*txt") into fusionout
  script:
  """
  bash ${repoDir}/process_scripts/alignment/starfusion.sh -p ${sampleid} -r ${index_path} -a ${fq1} -b ${fq2} -f
  """
}
process bamct {
  errorStrategy 'ignore'
  publishDir "$params.output/$caseid/rnaout", mode: 'copy'
  label 'profiling_qc'
  input:
  set caseid,sampleid,file(rbam),file(ridx) from ctbam
  output:
  file("${sampleid}.bamreadcount.txt.gz") into ctreads
  when:
  params.bamct == "detect"
  script:
  """
  export PATH=/project/shared/bicf_workflow_ref/seqprg/bam-readcount/bin/:$PATH
  bam-readcount -w 0 -q 0 -b 25 -f ${index_path}/genome.fa ${rbam} > ${sampleid}.bamreadcount.txt
  pigz ${sampleid}.bamreadcount.txt
  """
}
process alignqc {
  executor 'local'
  errorStrategy 'ignore'
  label 'profiling_qc'
  publishDir "$params.output/$caseid/rnaout", mode: 'copy'

  input:
  set caseid,sampleid,file(bam),file(hsout) from qcbam
  
  output:
  set file("${sampleid}_fastqc.zip"),file("${sampleid}_fastqc.html") into fastqc
  file("${sampleid}.sequence.stats.txt") into alignstats
  script:
  """
  bash ${repoDir}/process_scripts/alignment/bamqc.sh -p ${sampleid} -b ${bam} -n rna -e ${params.version}
  """
}
process geneabund {
  errorStrategy 'ignore'
  executor 'local'
  label 'geneabund'
  publishDir "$params.output/$caseid/rnaout", mode: 'copy'
  input:
  set caseid,sampleid, file(sbam) from abundbam
  output:
  file("${sampleid}.cts")  into counts
  file("${sampleid}_stringtie") into strcts
  file("${sampleid}.fpkm.txt") into fpkm
  file("${sampleid}.exonskip.answer.txt") into exonskip
  """
  bash ${repoDir}/process_scripts/genect_rnaseq/geneabundance.sh -s ${params.stranded} -g ${gtf_file} -p ${sampleid} -b ${sbam} -i ${ginfo} ${glist}
  bash ${repoDir}/process_scripts/genect_rnaseq/exonskipping.sh -g ${gtf_file} -p ${sampleid} -b ${sbam} -r ${index_path}
  """
}
process fb {
  queue '32GB'
  errorStrategy 'ignore'
  label 'variantcalling'
  publishDir "$params.output/$caseid/rnavcf", mode: 'copy'

  input:
  set caseid,$sampleid,file(gbam),file(gidx) from fbbam
  output:
  set caseid,file("${caseid}.fb*vcf.gz") into fbvcf
  script:
  """
  export biohpc=1
  bash ${repoDir}/process_scripts/variants/germline_vc.sh -r ${index_path} -p ${caseid} -a fb
  bash ${repoDir}/process_scripts/variants/uni_norm_annot.sh -g ${snpeff_vers} -r ${index_path} -p ${caseid}.fb -v ${caseid}.fb.vcf.gz 
  """
}
