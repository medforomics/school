#!/usr/bin/env nextflow

params.input = './analysis'
params.output = './analysis'

params.bams="$params.input/*.bam"
params.design="$params.input/design_tumor_only.txt"
params.genome="/project/shared/bicf_workflow_ref/human/grch38_cloud/dnaref"

index_path = file(params.genome)
design_file = file(params.design)
bams=file(params.bams)

params.cancer="detect"
params.callsvs="detect"
params.nuctype='dna'
params.projectid=''
snpeff_vers = 'GRCh38.86';

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

def oribam = []
def tarbam = []

new File(params.design).withReader { reader ->
    def hline = reader.readLine()
    def header = hline.split("\t")
    fidx = header.findIndexOf{it == 'FamilyID'};
    sidx = header.findIndexOf{it == 'SubjectID'};
    tidx = header.findIndexOf{it == 'SampleID'};
    oneidx = header.findIndexOf{it == 'BAM'};
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
	      oribam << tuple(row[fidx],row[tidx],fileMap.get(row[oneidx]))
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
   .into { fbbam; platbam; strelkabam; pindelbam; dellybam; svababam;}

gtxbam
   .groupTuple(by:0)		
   .set { gatkbam }

process pindel {
  errorStrategy 'ignore'
  queue '128GB,256GB,256GBv1'
  publishDir "$params.output/$subjid", mode: 'copy'
  input:
  set subjid,file(ssbam),file(ssidx) from pindelbam
  output:
  file("${subjid}.pindel_tandemdup.pass.vcf.gz") into tdvcf
  set subjid,file("${subjid}.pindel.vcf.gz") into pindelvcf
  file("${subjid}.pindel_sv.pass.vcf.gz") into pindelsv
  file("${subjid}.dna.genefusion.txt") into gf
  when:
  params.nuctype == "dna"
  script:
  """
  source /etc/profile.d/modules.sh
  module load samtools/gcc/1.8 snpeff/4.3q htslib/gcc/1.8 
  bash $baseDir/process_scripts/variants/svcalling.sh -r $index_path -b $ssbam -p $subjid -l ${index_path}/itd_genes.bed -a pindel
  perl $baseDir/process_scripts/variants/filter_pindel.pl -d ${subjid}.pindel_tandemdup.vcf.gz -s ${subjid}.pindel_sv.vcf.gz -i ${subjid}.pindel_indel.vcf.gz
  bgzip ${subjid}.pindel_indel.pass.vcf
  bgzip ${subjid}.pindel_tandemdup.pass.vcf
  mv ${subjid}.pindel_indel.pass.vcf.gz ${subjid}.pindel.vcf.gz
  grep '#CHROM' ${subjid}.pindel_sv.pass.vcf > ${subjid}.dna.genefusion.txt
  cat ${subjid}.pindel_sv.pass.vcf | \$SNPEFF_HOME/scripts/vcfEffOnePerLine.pl |java -jar \$SNPEFF_HOME/SnpSift.jar extractFields - CHROM POS END ANN[*].EFFECT ANN[*].GENE ANN[*].HGVS_C ANN[*].HGVS_P GEN[*] |grep -E 'CHROM|gene_fusion' |uniq >> ${subjid}.dna.genefusion.txt
  bgzip ${subjid}.pindel_sv.pass.vcf
  """
}
process delly {
  queue '32GB'
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid/${params.nuctype}_${params.projectid}", mode: 'copy'

  input:
  set subjid,file(ssbam),file(ssidx) from dellybam
  output:
  set subjid,file("${subjid}.delly.vcf.gz") into dellyvcf
  when:
  params.nuctype == "dna"
  script:				       
  """
  bash $baseDir/process_scripts/variants/svcalling.sh -r $index_path -b $ssbam -p $subjid -a delly
  """
}
process svaba {
  queue '32GB'
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid/${params.nuctype}_${params.projectid}", mode: 'copy'

  input:
  set subjid,file(ssbam),file(ssidx) from svababam
  output:
  set subjid,file("${subjid}.svaba.vcf.gz") into svabavcf
  set subjid,file("${subjid}.svaba.sv.vcf.gz") into svabasv
  when:
  params.nuctype == "dna"
  script:				       
  """
  bash $baseDir/process_scripts/variants/svcalling.sh -r $index_path -b $ssbam -p $subjid -a svaba
  """
}

process freebayes {
  queue '32GB'
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid/${params.nuctype}_${params.projectid}", mode: 'copy'

  input:
  set subjid,file(gbam),file(gidx) from fbbam
  output:
  set subjid,file("${subjid}.fb.vcf.gz") into fbvcf
  set subjid,file("${subjid}.fb.ori.vcf.gz") into fbori
  script:
  """
  bash $baseDir/process_scripts/variants/germline_vc.sh -r $index_path -p $subjid -a freebayes
  bash $baseDir/process_scripts/variants/uni_norm_annot.sh -g $snpeff_vers -r $index_path -p ${subjid}.fb -v ${subjid}.freebayes.vcf.gz 
  """
}
process mutect2 {
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
  bash $baseDir/process_scripts/variants/germline_vc.sh $ponopt -r $index_path -p $subjid -a mutect2
  bash $baseDir/process_scripts/variants/uni_norm_annot.sh -g $snpeff_vers -r $index_path -p ${subjid}.mutect -v ${subjid}.mutect.vcf.gz
  """
}
process strelka {
  queue '32GB'
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid/${params.nuctype}_${params.projectid}", mode: 'copy'

  input:
  set subjid,file(gbam),file(gidx) from strelkabam
  output:
  set subjid,file("${subjid}.strelka2.vcf.gz") into strelkavcf
  set subjid,file("${subjid}.strelka2.ori.vcf.gz") into strelkaori
  when:
  params.nuctype == "dna"
  script:
  """
  bash $baseDir/process_scripts/variants/germline_vc.sh -r $index_path -p $subjid -a strelka2
  bash $baseDir/process_scripts/variants/uni_norm_annot.sh -g $snpeff_vers -r $index_path -p ${subjid}.strelka2 -v ${subjid}.strelka2.vcf.gz
  """
}
process platypus {
  queue '32GB'
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid/${params.nuctype}_${params.projectid}", mode: 'copy'

  input:
  set subjid,file(gbam),file(gidx) from platbam
  output:
  set subjid,file("${subjid}.platypus.vcf.gz") into platvcf
  set subjid,file("${subjid}.platypus.ori.vcf.gz") into platori
  when:					       
  script:				       
  if (params.nuctype == "dna")
  """
  bash $baseDir/process_scripts/variants/germline_vc.sh -r $index_path -p $subjid -a platypus
  bash $baseDir/process_scripts/variants/uni_norm_annot.sh -g $snpeff_vers -r $index_path -p ${subjid}.platypus -v ${subjid}.platypus.vcf.gz
  """
  else
  """
  source /etc/profile.d/modules.sh
  module load samtools/1.6
  cp ${index_path}/union.header.vcf ${subjid}.platypus.vcf
  bgzip ${subjid}.platypus.vcf
  cp ${subjid}.platypus.vcf.gz ${subjid}.platypus.ori.vcf.gz
  """  
}

Channel
  .empty()
  .mix(fbvcf,platvcf,strelkavcf,gatkvcf,pindelvcf,svabavcf)
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
