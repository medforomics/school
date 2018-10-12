#!/usr/bin/env nextflow

params.input = './analysis'
params.output = './analysis'

params.bams="$params.input/*.bam"
params.design="$params.input/design.txt"

params.genome="/project/shared/bicf_workflow_ref/GRCh38"
params.targetpanel="$params.genome/clinseq_prj/UTSWV2.bed"
params.cancer="detect"
params.callsvs="detect"
params.nuctype='dna'
params.projectid=''

dbsnp_file="$params.genome/dbSnp.vcf.gz"
indel="$params.genome/GoldIndels.vcf.gz"
cosmic="$params.genome/cosmic.vcf.gz"
reffa=file("$params.genome/genome.fa")
capture=file(params.targetpanel)
design_file = file(params.design)
bams=file(params.bams)

dbsnp=file(dbsnp_file)
knownindel=file(indel)
index_path = file(params.genome)

snpeff_vers = 'GRCh38.86';

def fileMap = [:]

bams.each {
    final fileName = it.getFileName().toString()
    prefix = fileName.lastIndexOf('/')
    fileMap[fileName] = it
}

def oribam = []
def tarbam = []

new File(params.design).withReader { reader ->
    def hline = reader.readLine()
    def header = hline.split("\t")
    fidx = header.findIndexOf{it == 'FamilyID'};
    sidx = header.findIndexOf{it == 'SubjectID'};
    tidx = header.findIndexOf{it == 'SampleID'};
    oneidx = header.findIndexOf{it == 'BAM'};
    taridx = header.findIndexOf{it == 'FinalBAM'};
    if (sidx == -1) {
       sidx = tidx
    }
    if (fidx == -1) {
       fidx = sidx
    }
    while (line = reader.readLine()) {
    	   def row = line.split("\t")
	   if (fileMap.get(row[oneidx]) != null) {
	      oribam << tuple(row[fidx],row[tidx],fileMap.get(row[oneidx]))
	      tarbam << tuple(row[fidx],fileMap.get(row[taridx]))
	   }
	  
} 
}

if( ! oribam) { error "Didn't match any input files with entries in the design file" }
if( ! tarbam) { error "Didn't match any input files with entries in the design file" }

process indexoribams {
  errorStrategy 'ignore'
  input:
  set sid,tid,file(tumor) from oribam
  output:
  set sid,tid,file(tumor),file("${tumor}.bai") into svbam
  set sid,tid,file(tumor),file("${tumor}.bai") into cnvbam
  script:
  """
  bash $baseDir/process_scripts/alignment/indexbams.sh 
  """
}

process indexbams {
  input:
  set sid,file(tumor) from tarbam
  output:
  set sid,file(tumor),file("${tumor}.bai") into idxbam
  script:
  """
  bash $baseDir/process_scripts/alignment/indexbams.sh 
  """
}

idxbam
   .groupTuple(by:0)		
   .into { ssbam; sambam; hsbam; platbam; strelkabam }

process cnv {
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid/$pair_id", mode: 'copy'
  input:
  set subjid,pair_id,file(sbam),file(sidx) from cnvbam
  when:
  params.nuctype == "dna"
  output:
  file("${pair_id}.call.cns") into cns
  file("${pair_id}.cns") into cnsori
  file("${pair_id}.cnr") into cnr
  file("${pair_id}.answerplot*") into cnvansplot
  file("${pair_id}.*txt") into cnvtxt
  file("${pair_id}.cnv*pdf") into cnvpdf
  script:
  """
  bash $baseDir/process_scripts/variants/cnvkit.sh -u -c $capture -b $sbam -p $pair_id
  """
}
process pindel {
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid/$pair_id", mode: 'copy'
  input:
  set subjid,pair_id,file(ssbam),file(ssidx) from svbam
  output:
  file("${pair_id}.pindel_*.vcf.gz") into pindelvcf
  when:
  params.nuctype == "dna"
  script:
  """
  bash $baseDir/process_scripts/variants/pindel.sh -r ${index_path} -p ${pair_id}
  """
}

process mpileup {
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid/$params.nuctype", mode: 'copy'
  input:
  set subjid,file(gbam),file(gidx) from sambam
  
  output:
  set subjid,file("${subjid}.sam.vcf.gz") into samvcf
  set subjid,file("${subjid}.sam.ori.vcf.gz") into samori
  set subjid,file("${subjid}.sam.annot.vcf.gz") into samannot
  script:
  """
  bash $baseDir/process_scripts/variants/germline_vc.sh -r $index_path -p $subjid -a mpileup
  bash $baseDir/process_scripts/variants/norm_annot.sh -r $index_path -p ${subjid}.sam -v ${subjid}.sam.vcf.gz
  mv ${subjid}.sam.vcf.gz ${subjid}.sam.ori.vcf.gz
  mv ${subjid}.sam.norm.vcf.gz ${subjid}.sam.vcf.gz
  bash $baseDir/process_scripts/variants/annotvcf.sh -p ${subjid}.sam -r $index_path -v ${subjid}.sam.vcf.gz
  """
}
process hotspot {
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid/$params.nuctype", mode: 'copy'
  input:
  set subjid,file(gbam),file(gidx) from hsbam
  output:
  set subjid,file("${subjid}.hotspot.vcf.gz") into hsvcf
  set subjid,file("${subjid}.hotspot.ori.vcf.gz") into hsori
  set subjid,file("${subjid}.hotspot.annot.vcf.gz") into hsannot
  when:
  params.cancer == "detect"
  script:
  """
  bash $baseDir/process_scripts/variants/germline_vc.sh -r $index_path -p $subjid -a hotspot
  bash $baseDir/process_scripts/variants/norm_annot.sh -r $index_path -p ${subjid}.hotspot -v ${subjid}.hotspot.vcf.gz
  mv ${subjid}.hotspot.vcf.gz ${subjid}.hotspot.ori.vcf.gz
  mv ${subjid}.hotspot.norm.vcf.gz ${subjid}.hotspot.vcf.gz
  bash $baseDir/process_scripts/variants/annotvcf.sh -p ${subjid}.hotspot -r $index_path -v ${subjid}.hotspot.vcf.gz
  """
}
process speedseq {
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid/$params.nuctype", mode: 'copy'

  input:
  set subjid,file(gbam),file(gidx) from ssbam
  output:
  set subjid,file("${subjid}.ssvar.vcf.gz") into ssvcf
  set subjid,file("${subjid}.ssvar.ori.vcf.gz") into ssori
  set subjid,file("${subjid}.ssvar.annot.vcf.gz") into ssannot

  script:
  """
  bash $baseDir/process_scripts/variants/germline_vc.sh -r $index_path -p $subjid -a speedseq
  bash $baseDir/process_scripts/variants/norm_annot.sh -r $index_path -p ${subjid}.ssvar -v ${subjid}.ssvar.vcf.gz
  mv ${subjid}.ssvar.vcf.gz ${subjid}.ssvar.ori.vcf.gz
  mv ${subjid}.ssvar.norm.vcf.gz ${subjid}.ssvar.vcf.gz
  bash $baseDir/process_scripts/variants/annotvcf.sh -p ${subjid}.ssvar -r $index_path -v ${subjid}.ssvar.vcf.gz
  """
}

process strelka2 {
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid/$params.nuctype", mode: 'copy'

  input:
  set subjid,file(gbam),file(gidx) from strelkabam
  output:
  set subjid,file("${subjid}.strelka2.vcf.gz") into strelkavcf
  set subjid,file("${subjid}.strelka2.ori.vcf.gz") into strelkaori
  set subjid,file("${subjid}.strelka2.annot.vcf.gz") into strelkaannot
  script:
  if (params.nuctype == "dna")
  """
  bash $baseDir/process_scripts/variants/germline_vc.sh -r $index_path -p $subjid -a strelka2
  bash $baseDir/process_scripts/variants/norm_annot.sh -r $index_path -p ${subjid}.strelka2 -v ${subjid}.strelka2.vcf.gz
  mv ${subjid}.strelka2.vcf.gz ${subjid}.strelka2.ori.vcf.gz
  mv ${subjid}.strelka2.norm.vcf.gz ${subjid}.strelka2.vcf.gz
  bash $baseDir/process_scripts/variants/annotvcf.sh -p ${subjid}.strelka2 -r $index_path -v ${subjid}.strelka2.vcf.gz
  """
  else
  """
  source /etc/profile.d/modules.sh
  bash $baseDir/process_scripts/variants/germline_vc.sh -r $index_path -p $subjid -a gatk
  mv ${subjid}.gatk.vcf.gz ${subjid}.strelka2.vcf.gz
  bash $baseDir/process_scripts/variants/norm_annot.sh -r $index_path -p ${subjid}.strelka2 -v ${subjid}.strelka2.vcf.gz
  mv ${subjid}.strelka2.vcf.gz ${subjid}.strelka2.ori.vcf.gz
  mv ${subjid}.strelka2.norm.vcf.gz ${subjid}.strelka2.vcf.gz
  bash $baseDir/process_scripts/variants/annotvcf.sh -p ${subjid}.strelka2 -r $index_path -v ${subjid}.strelka2.vcf.gz
  """  
}

process platypus {
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid/$params.nuctype", mode: 'copy'

  input:
  set subjid,file(gbam),file(gidx) from platbam
  output:
  set subjid,file("${subjid}.platypus.vcf.gz") into platvcf
  set subjid,file("${subjid}.platypus.ori.vcf.gz") into platori
  set subjid,file("${subjid}.platypus.annot.vcf.gz") into platannot	
  when:					       
  script:				       
  if (params.nuctype == "dna")
  """
  bash $baseDir/process_scripts/variants/germline_vc.sh -r $index_path -p $subjid -a platypus
  bash $baseDir/process_scripts/variants/norm_annot.sh -r $index_path -p ${subjid}.platypus -v ${subjid}.platypus.vcf.gz
  mv ${subjid}.platypus.vcf.gz ${subjid}.platypus.ori.vcf.gz
  mv ${subjid}.platypus.norm.vcf.gz ${subjid}.platypus.vcf.gz
  bash $baseDir/process_scripts/variants/annotvcf.sh -p ${subjid}.platypus -r $index_path -v ${subjid}.platypus.vcf.gz
  """
  else
  """
  source /etc/profile.d/modules.sh
  module load samtools/1.6
  cp ${index_path}/union.header.vcf ${subjid}.platypus.vcf
  bgzip ${subjid}.platypus.vcf
  cp ${subjid}.platypus.vcf.gz ${subjid}.platypus.ori.vcf.gz
  cp ${subjid}.platypus.vcf.gz ${subjid}.platypus.annot.vcf.gz
  """  
}

Channel
  .empty()
  .mix(ssvcf,strelkavcf,samvcf,platvcf,hsvcf)
  .groupTuple(by:0)
  .set { vcflist}

process integrate {
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid", mode: 'copy'
  input:
  set subjid,file(vcf) from vcflist
  
  output:
  
  file("${subjid}${params.projectid}.*union.vcf.gz") into annotunionvcf
  file("${subjid}${params.projectid}.germline*vcf.gz") into annotvcf

  script:
  if (params.nuctype == "dna")
  """
  source /etc/profile.d/modules.sh
  bash $baseDir/process_scripts/variants/union.sh -r $index_path -p $subjid
  bash $baseDir/process_scripts/variants/annotvcf.sh -p $subjid -r $index_path -v ${subjid}.union.vcf.gz
  mv ${subjid}.union.vcf.gz ${subjid}${params.projectid}.dnaunion.vcf.gz
  mv ${subjid}.annot.vcf.gz ${subjid}${params.projectid}.germline.vcf.gz
  """
  else
  """
  source /etc/profile.d/modules.sh
  bash $baseDir/process_scripts/variants/union.sh -r $index_path -p $subjid
  bash $baseDir/process_scripts/variants/annotvcf.sh -p $subjid -r $index_path -v ${subjid}.union.vcf.gz
  mv ${subjid}.union.vcf.gz ${subjid}${params.projectid}.rnaunion.vcf.gz
  mv ${subjid}.annot.vcf.gz ${subjid}${params.projectid}.germline.rna.vcf.gz
  """
}
