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
fpalgo = ['fb', 'platypus']
svalgo = ['delly', 'svaba']
ssalgo = ['strelka2','shimmer']
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

process dtrim_align {
  queue '32GB,super'
  errorStrategy 'ignore'
  publishDir "$params.output/$caseid/$pair_id", mode: 'copy'
  input:
  set caseid,pair_id, file(fqs) from reads
  set caseid,tid,nid from ids
  output:
  set caseid,tid,nid,pair_id, file("${pair_id}.bam") into mdupbam
  set caseid,pair_id, file("${pair_id}.bam"),file("${pair_id}.bam.bai"),file("${pair_id}.trimreport.txt") into qcbam
  set caseid,pair_id, file("${pair_id}.bam") into virusalign
  set caseid,pair_id, file("${pair_id}.bam"),file("${pair_id}.bam.bai") into cnvbam
  set caseid,pair_id, file("${pair_id}.bam"),file("${pair_id}.bam.bai") into itdbam
  set caseid,tid,nid,file("${pair_id}.bam"), file("${pair_id}.bam.bai") into oribam
  script:
  """
  bash $baseDir/process_scripts/preproc_fastq/trimgalore.sh -f -p ${pair_id} ${fqs}
  bash $baseDir/process_scripts/alignment/dnaseqalign.sh -r $index_path -p $pair_id -x ${pair_id}.trim.R1.fastq.gz -y ${pair_id}.trim.R2.fastq.gz $alignopts
  """
}
process valign {
  executor 'local'
  errorStrategy 'ignore'
  publishDir "$params.output/$caseid/$pair_id", mode: 'copy'

  input:
  set caseid,pair_id, file(sbam) from virusalign
  output:
  file("${pair_id}.viral.seqstats.txt") into viralseqstats
  script:
  """
  bash $baseDir/process_scripts/alignment/virusalign.sh -b ${pair_id}.bam -p ${pair_id} -r $virus_index_path -f
  """
}

process markdups {
  errorStrategy 'ignore'
  queue '32GB,super'
  publishDir "$params.output/$caseid/$pair_id", mode: 'copy'

  input:
  set caseid,tid,nid,pair_id, file(sbam) from mdupbam
  output:
  set caseid,tid,nid,pair_id, file("${pair_id}.consensus.bam"),file("${pair_id}.consensus.bam.bai") into togatkbam
  set caseid,tid,nid,file("${pair_id}.consensus.bam"),file("${pair_id}.consensus.bam.bai") into consbam
  script:
  """
  bash $baseDir/process_scripts/alignment/markdups.sh -a $params.markdups -b $sbam -p $pair_id -r $index_path

  mv ${pair_id}.dedup.bam ${pair_id}.consensus.bam
  mv ${pair_id}.dedup.bam.bai ${pair_id}.consensus.bam.bai
  """
}

process dna_bamqc {
  errorStrategy 'ignore'
  publishDir "$params.output/$caseid/$pair_id", mode: 'copy'
  queue '128GB,256GB,256GBv1'
  input:
  set caseid,pair_id, file(gbam),file(idx),file(trimreport) from qcbam
  output:
  file("*fastqc*") into fastqc
  file("${pair_id}*txt") into dalignstats	
  script:
  """
  bash $baseDir/process_scripts/alignment/bamqc.sh -c $capturebed -n dna -r $index_path -b ${gbam} -p $pair_id -e ${params.version} 
  """
}

process cnv {
  executor 'local'
  errorStrategy 'ignore'
  publishDir "$params.output/$caseid/$pair_id", mode: 'copy'
  input:
  set caseid,pair_id,file(sbam),file(sidx) from cnvbam

  output:
  file("${pair_id}.call.cns") into cns
  file("${pair_id}.cns") into cnsori
  file("${pair_id}.cnr") into cnr
  file("${pair_id}.answerplot*") into cnvansplot
  file("${pair_id}.*txt") into cnvtxt
  file("${pair_id}.cnv*pdf") into cnvpdf
  when:
  skipCNV == false
  script:
  """
  bash $baseDir/process_scripts/variants/cnvkit.sh -r $index_path -b $sbam -p $pair_id -d $capturedir
  """
}

process itdseek {
  executor 'local'
  errorStrategy 'ignore'
  publishDir "$params.output/$caseid/$pair_id", mode: 'copy'
  input:
  set caseid,pair_id,file(sbam),file(sidx) from itdbam

  output:
  file("${pair_id}.itdseek_tandemdup.vcf.gz") into itdseekvcf

  script:
  """
  bash $baseDir/process_scripts/variants/svcalling.sh -b $sbam -r $index_path -p $pair_id -l ${index_path}/itd_genes.bed -a itdseek -g $params.snpeff_vers -f
  """
}

process gatkbam {
  queue '32GB,super'
  publishDir "$params.output/$caseid/$pair_id", mode: 'copy'

  input:
  set caseid,tid,nid,pair_id, file(sbam),file(idx) from togatkbam
  output:
  set caseid,tid,nid,file("${pair_id}.final.bam"),file("${pair_id}.final.bam.bai") into gtxbam
  script:
  """
  bash $baseDir/process_scripts/variants/gatkrunner.sh -a gatkbam -b $sbam -r $index_path -p $pair_id
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
  executor 'local'
  publishDir "$params.output/$caseid/dna_$params.seqrunid", mode: 'copy'
  errorStrategy 'ignore'
  input:
  set caseid,tid,nid,file(ssbam),file(ssidx) from msibam
  output:
  file("${caseid}*") into msiout
  script:
  if ( somatic[caseid] == true )
  """
  bash $baseDir/process_scripts/variants/msisensor.sh -r ${index_path} -p $caseid -b ${tid}.bam -n ${nid}.bam -c $capturebed
  """
  else
  """
  bash $baseDir/process_scripts/variants/msisensor.sh -r ${index_path} -p $caseid -b ${tid}.bam -c $capturebed
  """
}

process checkmates {
  executor 'local'
  publishDir "$params.output/$caseid/dna_$params.seqrunid", mode: 'copy'
  errorStrategy 'ignore'
  input:
  set caseid,tid,nid,file(bam),file(bidx) from checkbams
  output:
  file("${caseid}*") into checkmateout
  when: somatic[caseid] == true
  script:
  """
  bash $baseDir/process_scripts/variants/checkmate.sh -r ${index_path} -p ${pair_id} -c ${index_path}/NGSCheckMate.bed -f
  """
}

process pindel {
  errorStrategy 'ignore'
  queue '128GB,256GB,256GBv1'
  publishDir "$params.output/$caseid/dna_$params.seqrunid", mode: 'copy'
  input:
  set caseid,tid,nid,file(ssbam),file(ssidx) from pindelbam
  output:
  file("${caseid}.pindel_tandemdup.vcf.gz") into tdvcf
  set caseid,file("${caseid}.pindel.vcf.gz") into pindelvcf
  file("${caseid}.pindel.genefusion.txt") into pindelgf
  script:
  """
  bash $baseDir/process_scripts/variants/svcalling.sh -r $index_path -p $caseid -l ${index_path}/itd_genes.bed -a pindel -c ${index_path}/itd_genes.bed -g $params.snpeff_vers -f
  """
}

process sv {
  queue '32GB,super'
  errorStrategy 'ignore'
  publishDir "$params.output/$caseid/dna_$params.seqrunid", mode: 'copy'

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
  bash $baseDir/process_scripts/variants/svcalling.sh -r $index_path -x ${tid} -y ${nid} -b ${tid}.bam -n ${nid}.bam -p $caseid -a ${algo} -g $params.snpeff_vers -f 
  """
  else 
  """
  bash $baseDir/process_scripts/variants/svcalling.sh -r $index_path -b ${tid}.bam -p $caseid -a ${algo} -g $params.snpeff_vers -f
  """
}

process mutect {
  queue '128GB,256GB,256GBv1'
  errorStrategy 'ignore'
  publishDir "$params.output/$caseid/dna_$params.seqrunid", mode: 'copy'

  input:
  set caseid,tid,nid,file(ssbam),file(ssidx) from mutectbam
  output:
  set caseid,file("${caseid}.mutect.vcf.gz") into mutectvcf
  set caseid,file("${caseid}.mutect.ori.vcf.gz") into mutectori
  script:
  if ( somatic[caseid] == true ) 
  """
  bash $baseDir/process_scripts/variants/somatic_vc.sh $ponopt -r $index_path -p $caseid -x $tid -y $nid -t ${tid}.final.bam -n ${nid}.final.bam -b $capturebed -a mutect
  bash $baseDir/process_scripts/variants/uni_norm_annot.sh -g $params.snpeff_vers -r $index_path -p ${caseid}.mutect -v ${caseid}.mutect.vcf.gz
  """
  else
  """
  bash $baseDir/process_scripts/variants/germline_vc.sh $ponopt -r $index_path -p $caseid -b $capturebed -a mutect
  bash $baseDir/process_scripts/variants/uni_norm_annot.sh -g $params.snpeff_vers -r $index_path -p ${caseid}.mutect -v ${caseid}.mutect.vcf.gz
  """
}

process somvc {
  publishDir "$params.output/$caseid/dna_$params.seqrunid", mode: 'copy'
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
  bash $baseDir/process_scripts/variants/somatic_vc.sh -r $index_path -p $caseid -x $tid -y $nid -n ${nid}.consensus.bam -t ${tid}.consensus.bam -a ${algo} -b $capturebed
  bash $baseDir/process_scripts/variants/uni_norm_annot.sh -g $params.snpeff_vers -r $index_path -p ${caseid}.${algo} -v ${caseid}.${algo}.vcf.gz
  """
}

process germvc {
  queue '32GB,super'
  errorStrategy 'ignore'
  publishDir "$params.output/$caseid/dna_$params.seqrunid", mode: 'copy'
  input:
  set caseid,tid,nid,file(gbam),file(gidx) from germbam
  each algo from fpalgo
  output:
  set caseid,file("${caseid}.${algo}.vcf.gz") into germvcf
  set caseid,file("${caseid}.${algo}.ori.vcf.gz") into germori
  script:
  """
  bash $baseDir/process_scripts/variants/germline_vc.sh -r $index_path -p $caseid -a ${algo} -b $capturebed
  bash $baseDir/process_scripts/variants/uni_norm_annot.sh -g $params.snpeff_vers -r $index_path -p ${caseid}.${algo} -v ${caseid}.${algo}.vcf.gz 
  """
}

process germstrelka {
  queue '32GB,super'
  errorStrategy 'ignore'
  publishDir "$params.output/$caseid/dna_$params.seqrunid", mode: 'copy'

  input:
  set caseid,tid,nid,file(gbam),file(gidx) from strelkabam
  output:
  set caseid,file("${caseid}.strelka2.vcf.gz") into strelkavcf
  set caseid,file("${caseid}.strelka2.ori.vcf.gz") into strelkaori
  when: 
  somatic[caseid] == false
  script:
  """
  bash $baseDir/process_scripts/variants/germline_vc.sh -r $index_path -p $caseid -a strelka2 -b $capturebed
  bash $baseDir/process_scripts/variants/uni_norm_annot.sh -g $params.snpeff_vers -r $index_path -p ${caseid}.strelka2 -v ${caseid}.strelka2.vcf.gz 
  """
}

Channel
  .empty()
  .mix(mutectvcf,ssvcf,pindelvcf,germvcf,strelkavcf)
  .groupTuple(by:0)
  .set { vcflist}

process integrate {
  executor 'local'
  errorStrategy 'ignore'
  publishDir "$params.output/$caseid", mode: 'copy'
  input:
  set caseid,file(vcf) from vcflist
  output:
  file("${caseid}_${params.seqrunid}.dna.vcf.gz") into unionvcf
  script:
  """
  source /etc/profile.d/modules.sh
  module load htslib/gcc/1.8
  bash $baseDir/process_scripts/variants/union.sh -r $index_path -p $caseid
  cp ${caseid}.union.vcf.gz ${caseid}_${params.seqrunid}.dna.vcf.gz
  """
}
