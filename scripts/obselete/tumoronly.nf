#!/usr/bin/env nextflow

params.input = './analysis'
params.output = './analysis'

params.bams="$params.input/*.bam"
params.design="$params.input/design_tumor_only.txt"
params.genome="/project/shared/bicf_workflow_ref/human/grch38_cloud/dnaref"
params.capture="/project/shared/bicf_workflow_ref/human/grch38_cloud/panels/UTSW_V4_pancancer/targetpanel.bed"

index_path = file(params.genome)
design_file = file(params.design)
capturebed=file(params.capture)
bams=file(params.bams)

params.cancer="detect"
params.callsvs="detect"
params.nuctype='dna'
params.projectid=''
snpeff_vers = 'GRCh38.86';
fpsalgo = ['fb', 'strelka2', 'platypus']
svalgo = ['delly', 'svaba']
ponopt=''
if (params.pon) {
   ponopt="-q $params.pon"
}

def fileMap = [:]

bams.each {
    final fileName = it.getFileName().toString()
    prefix = fileName.lastIndexOf('/')
    fileMap[fileName] = it
}

def consbam =[]
def oribam = []
def tarbam = []

new File(params.design).withReader { reader ->
    def hline = reader.readLine()
    def header = hline.split("\t")
    fidx = header.findIndexOf{it == 'FamilyID'};
    sidx = header.findIndexOf{it == 'SubjectID'};
    tidx = header.findIndexOf{it == 'SampleID'};
    oneidx = header.findIndexOf{it == 'CBAM'};
    twoidx = header.findIndexOf{it == 'BAM'};
    gtidx = header.findIndexOf{it == 'GATKBAM'};
    if (sidx == -1) {
       sidx = tidx
    }
    if (fidx == -1) {
       fidx = sidx
    }
    while (line = reader.readLine()) {
    	   def row = line.split("\t")
	   if (fileMap.get(row[oneidx]) != null) {
	      consbam << tuple(row[fidx],row[tidx],fileMap.get(row[oneidx]))
	      oribam << tuple(row[fidx],row[tidx],fileMap.get(row[twoidx]))
	      tarbam << tuple(row[fidx],fileMap.get(row[gtidx]))
	   }
	  
} 
}

if( ! oribam) { error "Didn't match any input files with entries in the design file" }
if( ! tarbam) { error "Didn't match any input files with entries in the design file" }

process indexoribams {
  executor 'local'
  errorStrategy 'ignore'
  input:
  set sid,tid,file(tumor) from oribam
  output:
  set sid,file(tumor),file("${tumor}.bai") into idxbam
  script:
  """
  bash $baseDir/process_scripts/alignment/indexbams.sh 
  """
}

process indexconsbams {
  executor 'local'
  errorStrategy 'ignore'
  input:
  set sid,tid,file(tumor) from consbam
  output:
  set sid,file(tumor),file("${tumor}.bai") into cidxbam
  script:
  """
  bash $baseDir/process_scripts/alignment/indexbams.sh 
  """
}

process indextarbams {
  executor 'local'
  input:
  set sid,file(tumor) from tarbam
  output:
  set sid,file(tumor),file("${tumor}.bai") into gtxbam
  script:
  """
  bash $baseDir/process_scripts/alignment/indexbams.sh 
  """
}

idxbam
   .groupTuple(by:0)		
   .into { dellybam; svababam; msibam;}

cidxbam
   .groupTuple(by:0)		
   .into { fbbam; platbam; strelkabam; fpsbam; pindelbam;}

gtxbam
   .groupTuple(by:0)		
   .set { gatkbam }

process msi {
  executor 'local'
  publishDir "$params.output/$subjid/dna_$params.projectid", mode: 'copy'
  errorStrategy 'ignore'
  input:
  set subjid,file(ssbam),file(ssidx) from msibam
  output:
  file("${subjid}*") into msiout
  when:
  params.nuctype == "dna"
  script:
  """
  bash $baseDir/process_scripts/variants/msisensor.sh -r ${index_path} -p $subjid -b $ssbam -c $capturebed
  """
}

process pindel {
  errorStrategy 'ignore'
  queue '128GB,256GB,256GBv1'
  publishDir "$params.output/$subjid/${params.nuctype}_${params.projectid}", mode: 'copy'
  input:
  set subjid,file(ssbam),file(ssidx) from pindelbam
  output:
  file("${subjid}.pindel_tandemdup.vcf.gz") into tdvcf
  set subjid,file("${subjid}.pindel.vcf.gz") into pindelvcf
  file("${subjid}.pindel.genefusion.txt") into pindelgf
  when:
  params.nuctype == "dna"
  script:
  """
  source /etc/profile.d/modules.sh
  module load samtools/gcc/1.8 snpeff/4.3q htslib/gcc/1.8 
  bash $baseDir/process_scripts/variants/svcalling.sh -r $index_path -b $ssbam -p $subjid -l ${index_path}/itd_genes.bed -a pindel -f
  """
}

process sv {
  queue '32GB'
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid/${params.nuctype}_${params.projectid}", mode: 'copy'

  input:
  set subjid,file(ssbam),file(ssidx) from svababam
  each algo from svalgo
  output:
  set subjid,file("${subjid}.${algo}.vcf.gz") into svabavcf
  set subjid,file("${subjid}.${algo}.sv.vcf.gz") into svabasv
  file("${subjid}.${algo}.genefusion.txt") into svabagf
  when:
  params.nuctype == "dna"
  script:				       
  """
  bash $baseDir/process_scripts/variants/svcalling.sh -r $index_path -b $ssbam -p $subjid -a ${algo} -f
  """
}

process fps {
  queue '32GB'
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid/${params.nuctype}_${params.projectid}", mode: 'copy'

  input:
  set subjid,file(gbam),file(gidx) from fpsbam
  each algo from fpsalgo
  output:
  set subjid,file("${subjid}.${algo}.vcf.gz") into fpsvcf
  set subjid,file("${subjid}.${algo}.ori.vcf.gz") into fpsori
  script:
  """
  bash $baseDir/process_scripts/variants/germline_vc.sh -r $index_path -p $subjid -a ${algo}
  bash $baseDir/process_scripts/variants/uni_norm_annot.sh -g $snpeff_vers -r $index_path -p ${subjid}.${algo} -v ${subjid}.${algo}.vcf.gz 
  """
}
process mutect {
  queue '128GB,256GB,256GBv1'
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid/${params.nuctype}_${params.projectid}", mode: 'copy'

  input:
  set subjid,file(gbam),file(gidx) from gatkbam
  output:
  set subjid,file("${subjid}.mutect.vcf.gz") into gatkvcf
  set subjid,file("${subjid}.mutect.ori.vcf.gz") into gatkori
  script:
  when:
  params.nuctype == "dna"
  """
  bash $baseDir/process_scripts/variants/germline_vc.sh $ponopt -r $index_path -p $subjid -a mutect
  bash $baseDir/process_scripts/variants/uni_norm_annot.sh -g $snpeff_vers -r $index_path -p ${subjid}.mutect -v ${subjid}.mutect.vcf.gz
  """
}

Channel
  .empty()
  .mix(fpsvcf,gatkvcf,pindelvcf)
  .groupTuple(by:0)
  .set { vcflist}

process integrate {
  executor 'local'
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid", mode: 'copy'
  input:
  set subjid,file(vcf) from vcflist
  output:
  file("${subjid}_${params.projectid}*na.vcf.gz") into annotvcf
  script:
  if (params.nuctype == "dna")
  """
  source /etc/profile.d/modules.sh
  bash $baseDir/process_scripts/variants/union.sh -r $index_path -p $subjid
  mv ${subjid}.union.vcf.gz ${subjid}_${params.projectid}.dna.vcf.gz
  """
  else
  """
  source /etc/profile.d/modules.sh
  bash $baseDir/process_scripts/variants/union.sh -r $index_path -p $subjid
  mv ${subjid}.union.vcf.gz ${subjid}_${params.projectid}.rna.vcf.gz
  """
}