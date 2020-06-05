#!/usr/bin/env nextflow

params.input = './fastq'
params.output = './analysis'
snpeff_vers = 'GRCh38.86';
params.genome="/project/shared/bicf_workflow_ref/human/grch38_cloud/dnaref"
params.virus_genome="/project/shared/bicf_workflow_ref/human_virus_genome/clinlab_idt_genomes"

params.markdups='fgbio_umi'
params.version = 'v4'

params.caseid = 'test'
params.seqrunid = 'runtest'
params.vcfid = "${params.caseid}.${params.seqrunid}"
params.tumorid = 'Tumor'

somatic = false
fpsalgo = ['fb', 'strelka2', 'platypus']
svalgo = ['delly', 'svaba']
ssalgo = ['strelka2','shimmer']

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
    somatic = true
    fpsalgo = ['fb', 'platypus']
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

process dtrim_align {
  executor 'local'
  errorStrategy 'ignore'
  publishDir "$params.output/$params.caseid/$pair_id", mode: 'copy'
  input:
  set pair_id, file(fqs) from reads
  output:
  set pair_id, file("${pair_id}.bam") into mdupbam
  set pair_id, file("${pair_id}.bam"),file("${pair_id}.bam.bai"),file("${pair_id}.trimreport.txt") into qcbam
  set pair_id, file("${pair_id}.bam") into virusalign
  set pair_id, file("${pair_id}.bam"),file("${pair_id}.bam.bai") into cnvbam
  set pair_id, file("${pair_id}.bam"),file("${pair_id}.bam.bai") into itdbam
  set env(subjid),file("${pair_id}.bam"), file("${pair_id}.bam.bai") into oribam
  script:
  """
  subjid=${params.caseid}
  bash $baseDir/process_scripts/preproc_fastq/trimgalore.sh -f 1 -p ${pair_id} ${fqs}
  bash $baseDir/process_scripts/alignment/dnaseqalign.sh -r $index_path -p $pair_id -x ${pair_id}.trim.R1.fastq.gz -y ${pair_id}.trim.R2.fastq.gz $alignopts
  #bash $baseDir/process_scripts/cvc/dna_trim_align.sh -p ${pair_id} -r $index_path -f $alignopts $fqs
  """
}

process valign {
  executor 'local'
  errorStrategy 'ignore'
  publishDir "$params.output/$params.caseid/$pair_id", mode: 'copy'

  input:
  set pair_id, file(sbam) from virusalign
  output:
  file("${pair_id}.viral.seqstats.txt") into viralseqstats
  script:
  """
  bash $baseDir/process_scripts/alignment/virusalign.sh -b ${pair_id}.bam -p ${pair_id} -r $virus_index_path -f
  """
}

process markdups_consensus {
  errorStrategy 'ignore'
  queue '32GB,super'
  publishDir "$params.output/$params.caseid/$pair_id", mode: 'copy'

  input:
  set  pair_id, file(sbam) from mdupbam
  output:
  set pair_id, file("${pair_id}.consensus.bam"),file("${pair_id}.consensus.bam.bai") into togatkbam
  set env(subjid),file("${pair_id}.consensus.bam"),file("${pair_id}.consensus.bam.bai") into consbam
  script:
  """
  subjid=${params.caseid}
  bash $baseDir/process_scripts/alignment/markdups.sh -a $params.markdups -b $sbam -p $pair_id -r $index_path

  mv ${pair_id}.dedup.bam ${pair_id}.consensus.bam
  mv ${pair_id}.dedup.bam.bai ${pair_id}.consensus.bam.bai
  """
}

process qc_gbam           {
  errorStrategy 'ignore'
  publishDir "$params.output/$params.caseid/$pair_id", mode: 'copy'
  queue '128GB,256GB,256GBv1'
  input:
  set pair_id, file(gbam),file(idx),file(trimreport) from qcbam
  output:
  file("*fastqc*") into fastqc
  file("${pair_id}*txt") into dalignstats	
  script:
  """
  bash $baseDir/process_scripts/alignment/bamqc.sh -c $capturebed -n dna -r $index_path -b ${gbam} -p $pair_id  
  perl $baseDir/scripts/sequenceqc_dna.pl -r ${index_path} -e ${params.version} *.genomecov.txt
  """
}

process cnv {
  executor 'local'
  errorStrategy 'ignore'
  publishDir "$params.output/$params.caseid/$pair_id", mode: 'copy'
  input:
  set pair_id,file(sbam),file(sidx) from cnvbam

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
  publishDir "$params.output/$params.caseid/$pair_id", mode: 'copy'
  input:
  set pair_id,file(sbam),file(sidx) from itdbam

  output:
  file("${pair_id}.itdseek_tandemdup.vcf.gz") into itdseekvcf

  script:
  """
  bash $baseDir/process_scripts/variants/svcalling.sh -b $sbam -r $index_path -p $pair_id -l ${index_path}/itd_genes.bed -a itdseek -f
  """
}

process gatkbam {
  queue '32GB,super'
  publishDir "$params.output/$params.caseid/$pair_id", mode: 'copy'

  input:
  set pair_id, file(sbam),file(idx) from togatkbam
  output:
  set env(subjid),file("${pair_id}.final.bam"),file("${pair_id}.final.bam.bai") into gtxbam
  script:
  """
  subjid=${params.caseid}
  bash $baseDir/process_scripts/variants/gatkrunner.sh -a gatkbam -b $sbam -r $index_path -p $pair_id
  """
}

oribam
   .groupTuple(by:0)
   .into { svbam; msibam;}

consbam
   .groupTuple(by:0)
   .into { checkbams; sombam; germbam; pindelbam;}

gtxbam
   .groupTuple(by:0)		
   .set { mutectbam }

process msi {
  executor 'local'
  publishDir "$params.output/$params.caseid/dna_$params.seqrunid", mode: 'copy'
  errorStrategy 'ignore'
  input:
  set caseid,file(ssbam),file(ssidx) from msibam
  output:
  file("${caseid}*") into msiout
  script:
  if ( somatic == true )
  """
  bash $baseDir/process_scripts/variants/msisensor.sh -r ${index_path} -p $params.caseid -b ${params.tumorid}.bam -n ${params.normalid}.bam -c $capturebed
  """
  else
  """
  bash $baseDir/process_scripts/variants/msisensor.sh -r ${index_path} -p $params.caseid -b ${params.tumorid}.bam -c $capturebed
  """
}

process checkmates {
  executor 'local'
  publishDir "$params.output/$params.caseid/dna_$params.seqrunid", mode: 'copy'
  errorStrategy 'ignore'
  input:
  set caseid,file(bam),file(bidx) from checkbams
  file(conf) from ncmconf
  output:
  file("${caseid}*") into checkmateout
  when: somatic == true
  script:
  """
  source /etc/profile.d/modules.sh
  module load python/2.7.x-anaconda git/v2.5.3
  python /project/shared/bicf_workflow_ref/seqprg/NGSCheckMate/ncm.py -B -d ./ -bed ${index_path}/NGSCheckMate.bed -O ./ -N ${caseid}
  perl $baseDir/scripts/sequenceqc_somatic.pl -r ${index_path} -i ${caseid}_all.txt -o ${caseid}_${params.seqrunid}.sequence.stats.txt
  """
}

process pindel {
  errorStrategy 'ignore'
  queue '128GB,256GB,256GBv1'
  publishDir "$params.output/$params.caseid/dna_$params.seqrunid", mode: 'copy'
  input:
  set caseid,file(ssbam),file(ssidx) from pindelbam
  output:
  file("${caseid}.pindel_tandemdup.vcf.gz") into tdvcf
  set caseid,file("${caseid}.pindel.vcf.gz") into pindelvcf
  file("${caseid}.pindel.genefusion.txt") into pindelgf
  script:
  """
  bash $baseDir/process_scripts/variants/svcalling.sh -r $index_path -p $params.caseid -l ${index_path}/itd_genes.bed -a pindel -f
  """
}

process sv {
  queue '32GB,super'
  errorStrategy 'ignore'
  publishDir "$params.output/$params.caseid/dna_$params.seqrunid", mode: 'copy'

  input:
  set caseid,file(ssbam),file(ssidx) from svbam
  each algo from svalgo
  output:
  set caseid,file("${caseid}.${algo}.vcf.gz") into svvcf
  set caseid,file("${caseid}.${algo}.sv.vcf.gz") optional true into svsv
  file("${caseid}.${algo}.genefusion.txt") into svgf

  script:				       
  if ( somatic == true ) 
  """
  bash $baseDir/process_scripts/variants/svcalling.sh -r $index_path -x ${params.tumorid} -y ${params.normalid} -b ${params.tumorid}.bam -n ${params.normalid}.bam -p $params.caseid -a ${algo} -f 
  """
  else 
  """
  bash $baseDir/process_scripts/variants/svcalling.sh -r $index_path -b ${params.tumorid}.bam -p $params.caseid -a ${algo} -f
  """
}

process mutect {
  queue '128GB,256GB,256GBv1'
  errorStrategy 'ignore'
  publishDir "$params.output/$params.caseid/dna_$params.seqrunid", mode: 'copy'

  input:
  set caseid,file(ssbam),file(ssidx) from svbam from mutectbam
  output:
  set caseid,file("${caseid}.mutect.vcf.gz") into mutectvcf
  set caseid,file("${caseid}.mutect.ori.vcf.gz") into mutectori
  script:
  if ( somatic == true ) 
  """
  bash $baseDir/process_scripts/variants/somatic_vc.sh $ponopt -r $index_path -p $params.caseid -x $params.tumorid -y $params.normalid -t ${params.tumorid}.final.bam -n ${params.normalid}.final.bam -a mutect
  bash $baseDir/process_scripts/variants/uni_norm_annot.sh -g $snpeff_vers -r $index_path -p ${caseid}.mutect -v ${caseid}.mutect.vcf.gz
  """
  else
  """
  bash $baseDir/process_scripts/variants/germline_vc.sh $ponopt -r $index_path -p $params.caseid -a mutect
  bash $baseDir/process_scripts/variants/uni_norm_annot.sh -g $snpeff_vers -r $index_path -p ${caseid}.mutect -v ${caseid}.mutect.vcf.gz
  """
}

process somvc {
  queue '32GB,super'
  errorStrategy 'ignore'
  publishDir "$params.output/$params.caseid/dna_$params.seqrunid", mode: 'copy'

  input:
  set caseid,file(ssbam),file(ssidx) from svbam from sombam
  each algo from ssalgo
  output:
  set caseid,file("${caseid}.${algo}.vcf.gz") into ssvcf
  set caseid,file("${caseid}.${algo}.ori.vcf.gz") into ssori
  when:
  somatic == true
  script:
  """
  bash $baseDir/process_scripts/variants/somatic_vc.sh -r $index_path -p $params.caseid -x $params.tumorid -y $params.normalid -n ${params.normalid}.consensus.bam -t ${params.tumorid}.consensus.bam -a ${algo} -b $capturebed
  bash $baseDir/process_scripts/variants/uni_norm_annot.sh -g $snpeff_vers -r $index_path -p ${caseid}.${algo} -v ${caseid}.${algo}.vcf.gz
  """
}

process germvc {
  queue '32GB,super'
  errorStrategy 'ignore'
  publishDir "$params.output/$params.caseid/dna_$params.seqrunid", mode: 'copy'


  input:
  set caseid,file(gbam),file(gidx) from germbam
  each algo from fpsalgo
  output:
  set caseid,file("${caseid}.${algo}.vcf.gz") into germvcf
  set caseid,file("${caseid}.${algo}.ori.vcf.gz") into germori
  script:
  """
  bash $baseDir/process_scripts/variants/germline_vc.sh -r $index_path -p $params.caseid -a ${algo} -b $capturebed
  bash $baseDir/process_scripts/variants/uni_norm_annot.sh -g $snpeff_vers -r $index_path -p ${caseid}.${algo} -v ${caseid}.${algo}.vcf.gz 
  """
}

Channel
  .empty()
  .mix(mutectvcf,ssvcf,pindelvcf,germvcf)
  .groupTuple(by:0)
  .set { vcflist}

process integrate {
  executor 'local'
  errorStrategy 'ignore'
  publishDir "$params.output/$params.caseid", mode: 'copy'
  input:
  set caseid,file(vcf) from vcflist
  output:
  file("${caseid}_${params.seqrunid}.dna.vcf.gz") into unionvcf
  script:
  """
  source /etc/profile.d/modules.sh
  module load htslib/gcc/1.8
  bash $baseDir/process_scripts/variants/union.sh -r $index_path -p $params.caseid
  cp ${caseid}.union.vcf.gz ${caseid}_${params.seqrunid}.dna.vcf.gz
  """
}
