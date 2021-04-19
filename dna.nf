#!/usr/bin/env nextflow

params.input = './fastq'
params.output = './analysis'
params.snpeff_vers = 'GRCh38.86';
params.genome="/project/shared/bicf_workflow_ref/human/grch38_cloud/dnaref"
params.virus_genome="/project/shared/bicf_workflow_ref/human_virus_genome/clinlab_idt_genomes"

params.markdups='fgbio_umi'
params.version = 'v4'
params.seqrunid = 'runtest'

somatic = false
fpalgo = ['fb']
ssalgo = ['strelka2']
svalgo = ['delly', 'svaba']

ncmconf = file("$params.genome/ncm.conf")
reffa=file("$params.genome/genome.fa")
dbsnp="$params.genome/dbSnp.vcf.gz"
indel="$params.genome/GoldIndels.vcf.gz"

dbsnp=file(dbsnp)
knownindel=file(indel)
index_path=file(params.genome)
capturebed = file("$params.capture")
capturedir = file("$params.capturedir")
virus_index_path=file(params.virus_genome)

skipCNV = false
if(capturedir.isEmpty()) {
  skipCNV = true
}
alignopts = ''
if (params.markdups == 'fgbio_umi') {
   alignopts='-u'
}
repoDir=workflow.projectDir
if (params.repoDir) {
   repoDir=params.repoDir
}
ponopt=''
if (params.pon) {
   ponopt="-q $params.pon"
}
def somatic = [:]
def ids = []
def reads = []

if (params.caseid) {
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
   oneidx = header.findIndexOf{it == 'FqR1'};
   twoidx = header.findIndexOf{it == 'FqR2'};
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

process dtrim {
  queue '32GB,super'
  label 'trim'
  errorStrategy 'ignore'
  publishDir "$params.output/$caseid/dnaout", mode: 'copy'
  input:
  set caseid,sampleid, file(fqs) from reads
  set caseid,tid,nid from ids
  output:
  set caseid,tid,nid,sampleid,file("${sampleid}.trim.R1.fastq.gz"),file("${sampleid}.trim.R2.fastq.gz"),file("${sampleid}.trimreport.txt") into treads
  script:
  """
  bash ${repoDir}/process_scripts/preproc_fastq/trimgalore.sh -f -p ${sampleid} ${fqs}
  """
}
process dalign {
  queue '32GB,super'
  label 'dnaalign'
  errorStrategy 'ignore'
  publishDir "$params.output/$caseid/dnaout", mode: 'copy'
  input:
  set caseid,tid,nid,sampleid,file(fq1),file(fq2),file(trimout) from treads
  output:
  set caseid,tid,nid,sampleid, file("${sampleid}.bam") into mdupbam
  set caseid,sampleid, file("${sampleid}.bam"),file("${sampleid}.bam.bai"),file(trimout) into qcbam
  set caseid,sampleid,file("${sampleid}.bam") into virusalign
  set caseid,sampleid,file("${sampleid}.bam"),file("${sampleid}.bam.bai") into cnvbam
  set caseid,sampleid,file("${sampleid}.bam"),file("${sampleid}.bam.bai") into itdbam
  set caseid,tid,nid,file("${sampleid}.bam"), file("${sampleid}.bam.bai") into oribam
  script:
  """
  bash ${repoDir}/process_scripts/alignment/dnaseqalign.sh -r $index_path -p $sampleid -x ${fq1} -y ${fq2} $alignopts
  """
}
process valign {
  queue '32GB,super'
  label 'dnaalign'
  errorStrategy 'ignore'
  publishDir "$params.output/$caseid/dnaout", mode: 'copy'

  input:
  set caseid,sampleid, file(sbam) from virusalign
  output:
  file("${sampleid}.viral.seqstats.txt") into viralseqstats
  script:
  """
  bash ${repoDir}/process_scripts/alignment/virusalign.sh -b ${sampleid}.bam -p ${sampleid} -r $virus_index_path -f
  """
}

process markdups {
  errorStrategy 'ignore'
  label 'dnaalign'
  queue '32GB,super'
  publishDir "$params.output/$caseid/dnaout", mode: 'copy'

  input:
  set caseid,tid,nid,sampleid, file(sbam) from mdupbam
  output:
  set caseid,tid,nid,sampleid, file("${sampleid}.consensus.bam"),file("${sampleid}.consensus.bam.bai") into togatkbam
  set caseid,tid,nid,file("${sampleid}.consensus.bam"),file("${sampleid}.consensus.bam.bai") into consbam
  script:
  """
  bash ${repoDir}/process_scripts/alignment/markdups.sh -a $params.markdups -b $sbam -p $sampleid -r $index_path

  mv ${sampleid}.dedup.bam ${sampleid}.consensus.bam
  mv ${sampleid}.dedup.bam.bai ${sampleid}.consensus.bam.bai
  """
}

process dna_bamqc {
  errorStrategy 'ignore'
  label 'profiling_qc'
  publishDir "$params.output/$caseid/dnaout", mode: 'copy'
  queue '128GB,256GB,256GBv1'
  input:
  set caseid,sampleid, file(gbam),file(idx),file(trimreport) from qcbam
  output:
  file("*fastqc*") into fastqc
  file("${sampleid}*txt") into dalignstats	
  script:
  """
  bash ${repoDir}/process_scripts/alignment/bamqc.sh -c $capturebed -n dna -r $index_path -b ${gbam} -p $sampleid -e ${params.version} 
  """
}

process cnv {
  queue '32GB,super'
  label 'structuralvariant'
  errorStrategy 'ignore'
  publishDir "$params.output/$caseid/dnacallset", mode: 'copy'
  input:
  set caseid,sampleid,file(sbam),file(sidx) from cnvbam

  output:
  file("${sampleid}.call.cns") into cns
  file("${sampleid}.cns") into cnsori
  file("${sampleid}.cnr") into cnr
  file("${sampleid}.answerplot*") into cnvansplot
  file("${sampleid}.*txt") into cnvtxt
  file("${sampleid}.cnv*pdf") into cnvpdf
  when:
  skipCNV == false
  script:
  """
  bash ${repoDir}/process_scripts/variants/cnvkit.sh -r $index_path -b $sbam -p $sampleid -d $capturedir
  """
}

process itdseek {
  queue '32GB,super'
  label 'structuralvariant'
  errorStrategy 'ignore'
  publishDir "$params.output/$caseid/dnacallset", mode: 'copy'
  input:
  set caseid,sampleid,file(sbam),file(sidx) from itdbam

  output:
  file("${sampleid}.itdseek_tandemdup.vcf.gz") into itdseekvcf

  script:
  """
  bash ${repoDir}/process_scripts/variants/svcalling.sh -b $sbam -r $index_path -p $sampleid -l ${index_path}/itd_genes.bed -a itdseek -g $params.snpeff_vers -f
  """
}

process gatkbam {
  queue '32GB,super'
  label 'variantcalling'
  publishDir "$params.output/$caseid/dnaout", mode: 'copy'

  input:
  set caseid,tid,nid,sampleid, file(sbam),file(idx) from togatkbam
  output:
  set caseid,tid,nid,file("${sampleid}.final.bam"),file("${sampleid}.final.bam.bai") into gtxbam
  script:
  """
  bash ${repoDir}/process_scripts/variants/gatkrunner.sh -a gatkbam -b $sbam -r $index_path -p $sampleid
  """
}

oribam
   .groupTuple(by:[0,1,2])
   .into { svbam; msibam; }

consbam
   .groupTuple(by:[0,1,2])
   .into { checkbams; sombam; germbam; pindelbam; strelkabam; }

gtxbam
   .groupTuple(by:[0,1,2])		
   .set { mutectbam }

process msi {
  queue '32GB,super'
  label 'profiling_qc'
  publishDir "$params.output/$caseid/dnacallset", mode: 'copy'
  errorStrategy 'ignore'
  input:
  set caseid,tid,nid,file(ssbam),file(ssidx) from msibam
  output:
  file("${caseid}*") into msiout
  script:
  if ( somatic[caseid] == true )
  """
  bash ${repoDir}/process_scripts/variants/msisensor.sh -r ${index_path} -p $caseid -b ${tid}.bam -n ${nid}.bam -c $capturebed
  """
  else
  """
  bash ${repoDir}/process_scripts/variants/msisensor.sh -r ${index_path} -p $caseid -b ${tid}.bam -c $capturebed
  """
}

process checkmates {
   queue '32GB,super'
  label 'profiling_qc'
  publishDir "$params.output/$caseid/dnacallset", mode: 'copy'
  errorStrategy 'ignore'
  input:
  set caseid,tid,nid,file(bam),file(bidx) from checkbams
  output:
  file("${caseid}*") into checkmateout
  when: somatic[caseid] == true
  script:
  """
  bash ${repoDir}/process_scripts/variants/checkmate.sh -r ${index_path} -p ${caseid} -c ${index_path}/NGSCheckMate.bed -f
  """
}

process pindel {
  errorStrategy 'ignore'
  label 'structuralvariant'
  queue '128GB,256GB,256GBv1'
  publishDir "$params.output/$caseid/dnacallset", mode: 'copy'
  input:
  set caseid,tid,nid,file(ssbam),file(ssidx) from pindelbam
  output:
  file("${caseid}.pindel_tandemdup.vcf.gz") into tdvcf
  set caseid,file("${caseid}.pindel.vcf.gz") into pindelvcf
  file("${caseid}.pindel.genefusion.txt") into pindelgf
  script:
  """
  bash ${repoDir}/process_scripts/variants/svcalling.sh -r $index_path -p $caseid -l ${index_path}/itd_genes.bed -a pindel -c ${index_path}/cancer_gene_list.bed -g $params.snpeff_vers -f
  """
}

process sv {
  queue '32GB,super'
  label 'structuralvariant'
  errorStrategy 'ignore'
  publishDir "$params.output/$caseid/dnacallset", mode: 'copy'

  input:
  set caseid,tid,nid,file(ssbam),file(ssidx) from svbam
  each algo from svalgo
  output:
  set caseid,file("${caseid}.${algo}.vcf.gz") into svvcf
  set caseid,file("${caseid}.${algo}.sv.vcf.gz") optional true into svsv
  file("${caseid}.${algo}.genefusion.txt") into svgf

  script:				       
  if ( somatic[caseid] == true ) 
  """
  bash ${repoDir}/process_scripts/variants/svcalling.sh -r $index_path -x ${tid} -y ${nid} -b ${tid}.bam -n ${nid}.bam -p $caseid -a ${algo} -g $params.snpeff_vers -f 
  """
  else 
  """
  bash ${repoDir}/process_scripts/variants/svcalling.sh -r $index_path -b ${tid}.bam -p $caseid -a ${algo} -g $params.snpeff_vers -f
  """
}

process mutect {
  queue '128GB,256GB,256GBv1'
  label 'variantcalling'
  errorStrategy 'ignore'
  publishDir "$params.output/$caseid/dnacallset", mode: 'copy'

  input:
  set caseid,tid,nid,file(ssbam),file(ssidx) from mutectbam
  output:
  set caseid,file("${caseid}.mutect.vcf.gz") into mutectvcf
  set caseid,file("${caseid}.mutect.ori.vcf.gz") into mutectori
  script:
  if ( somatic[caseid] == true ) 
  """
  bash ${repoDir}/process_scripts/variants/somatic_vc.sh $ponopt -r $index_path -p $caseid -x $tid -y $nid -t ${tid}.final.bam -n ${nid}.final.bam -b $capturebed -a mutect
  bash ${repoDir}/process_scripts/variants/uni_norm_annot.sh -g $params.snpeff_vers -r $index_path -p ${caseid}.mutect -v ${caseid}.mutect.vcf.gz
  """
  else
  """
  bash ${repoDir}/process_scripts/variants/germline_vc.sh $ponopt -r $index_path -p $caseid -b $capturebed -a mutect
  bash ${repoDir}/process_scripts/variants/uni_norm_annot.sh -g $params.snpeff_vers -r $index_path -p ${caseid}.mutect -v ${caseid}.mutect.vcf.gz
  """
}

process somvc {
  publishDir "$params.output/$caseid/dnacallset", mode: 'copy'
  label 'variantcalling'
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
  maxErrors 20
  queue '32GB,super'

  input:
  set caseid,tid,nid,file(ssbam),file(ssidx) from sombam
  each algo from ssalgo
  output:
  set caseid,file("${caseid}.${algo}.vcf.gz") into ssvcf
  set caseid,file("${caseid}.${algo}.ori.vcf.gz") into ssori
  when:
  somatic[caseid] == true
  script:
  """
  bash ${repoDir}/process_scripts/variants/somatic_vc.sh -r $index_path -p $caseid -x $tid -y $nid -n ${nid}.consensus.bam -t ${tid}.consensus.bam -a ${algo} -b $capturebed
  bash ${repoDir}/process_scripts/variants/uni_norm_annot.sh -g $params.snpeff_vers -r $index_path -p ${caseid}.${algo} -v ${caseid}.${algo}.vcf.gz
  """
}

process germvc {
  queue '32GB,super'
  label 'variantcalling'
  errorStrategy 'ignore'
  publishDir "$params.output/$caseid/dnacallset", mode: 'copy'
  input:
  set caseid,tid,nid,file(gbam),file(gidx) from germbam
  each algo from fpalgo
  output:
  set caseid,file("${caseid}.${algo}.vcf.gz") into germvcf
  set caseid,file("${caseid}.${algo}.ori.vcf.gz") into germori
  script:
  """
  bash ${repoDir}/process_scripts/variants/germline_vc.sh -r $index_path -p $caseid -a ${algo} -b $capturebed
  bash ${repoDir}/process_scripts/variants/uni_norm_annot.sh -g $params.snpeff_vers -r $index_path -p ${caseid}.${algo} -v ${caseid}.${algo}.vcf.gz 
  """
}

process germstrelka {
  queue '32GB,super'
  label 'variantcalling'
  errorStrategy 'ignore'
  publishDir "$params.output/$caseid/dnacallset", mode: 'copy'

  input:
  set caseid,tid,nid,file(gbam),file(gidx) from strelkabam
  output:
  set caseid,file("${caseid}.strelka2.vcf.gz") into strelkavcf
  set caseid,file("${caseid}.strelka2.ori.vcf.gz") into strelkaori
  when: 
  somatic[caseid] == false
  script:
  """
  bash ${repoDir}/process_scripts/variants/germline_vc.sh -r $index_path -p $caseid -a strelka2 -b $capturebed
  bash ${repoDir}/process_scripts/variants/uni_norm_annot.sh -g $params.snpeff_vers -r $index_path -p ${caseid}.strelka2 -v ${caseid}.strelka2.vcf.gz 
  """
}

Channel
  .empty()
  .mix(mutectvcf,ssvcf,pindelvcf,germvcf,strelkavcf)
  .groupTuple(by:0)
  .set { vcflist}

process integrate {
  executor 'local'
  label 'variantcalling'
  errorStrategy 'ignore'
  publishDir "$params.output/$caseid/dnavcf", mode: 'copy'
  input:
  set caseid,file(vcf) from vcflist
  output:
  file("${caseid}.union.vcf.gz") into unionvcf
  script:
  """
  bash ${repoDir}/process_scripts/variants/union.sh -r $index_path -p $caseid
  #cp ${caseid}.union.vcf.gz ${caseid}.dna.vcf.gz
  """
}
