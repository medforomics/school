#!/usr/bin/env nextflow

// Default parameter values to run tests
params.input = './'
params.output = './'

params.bams="$params.input/*.bam"
params.design="$params.input/design.txt"
params.callsvs="skip"
params.genome="/project/shared/bicf_workflow_ref/human/grch38_cloud/dnaref"
params.capture="/project/shared/bicf_workflow_ref/human/grch38_cloud/panels/UTSW_V4_pancancer/targetpanel.bed"
params.projectid=''

design_file = file(params.design)
bams=file(params.bams)
capturebed=file(params.capture)
index_path = file(params.genome)
ncmconf = file("$params.genome/ncm.conf")

fpsalgo = ['fb', 'platypus']
svalgo = ['delly', 'svaba']
ssalgo = ['strelka2','shimmer']

ponopt=''
if (params.pon) {
   ponopt="-q $params.pon"
}

snpeff_vers = 'GRCh38.86';

def fileMap = [:]

bams.each {
    final fileName = it.getFileName().toString()
    prefix = fileName.lastIndexOf('/')
    fileMap[fileName] = it
}
def oribam = []
def tarbam = []
def consbam = []
def nameMap = [:]

new File(params.design).withReader { reader ->
    def hline = reader.readLine()
    def header = hline.split("\t")
    pidx = header.findIndexOf{it == 'PairID'};
    vidx = header.findIndexOf{it == 'VcfID'};
    tidx = header.findIndexOf{it == 'TumorID'};
    nidx = header.findIndexOf{it == 'NormalID'};
    oneidx = header.findIndexOf{it == 'TumorBAM'};
    twoidx = header.findIndexOf{it == 'NormalBAM'};
    ctidx = header.findIndexOf{it == 'TumorCBAM'};
    cnidx = header.findIndexOf{it == 'NormalCBAM'};
    totidx = header.findIndexOf{it == 'TumorGATKBAM'};
    notidx = header.findIndexOf{it == 'NormalGATKBAM'};

    if (twoidx == -1) {
       twoidx = oneidx
       }      
    while (line = reader.readLine()) {
    	   def row = line.split("\t")
	   if (fileMap.get(row[oneidx]) != null) {
	      nameMap[row[pidx]] = row[vidx]
	      oribam << tuple(row[pidx],row[tidx],row[nidx],fileMap.get(row[oneidx]),fileMap.get(row[twoidx]))
	      consbam << tuple(row[pidx],row[tidx],row[nidx],fileMap.get(row[ctidx]),fileMap.get(row[cnidx]))
	      tarbam << tuple(row[pidx],row[tidx],row[nidx],fileMap.get(row[totidx]),fileMap.get(row[notidx]))
	   }
} 
}

if( ! oribam) { error "Didn't match any input files with entries in the design file" }
if( ! tarbam) { error "Didn't match any input files with entries in the design file" }

process indexoribams {
  executor 'local'
  errorStrategy 'ignore'
  input:
  set pid,tid,nid,file(tumor),file(normal) from oribam
  output:
  set pid,tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into idxbam
  script:
  """
  bash $baseDir/process_scripts/alignment/indexbams.sh 
  """
}

process indexconsbams {
  executor 'local'
  errorStrategy 'ignore'
  input:
  set pid,tid,nid,file(tumor),file(normal) from consbam
  output:
  set pid,tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into cidxbam
  script:
  """
  bash $baseDir/process_scripts/alignment/indexbams.sh 
  """
}

process indextarbams {
  executor 'local'
  errorStrategy 'ignore'
  input:
  set pid,tid,nid,file(tumor),file(normal) from tarbam
  output:
  set pid,tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into mutectbam
  set pid,tid,nid into pairnames
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
   .into { fbbam; platbam; strelkabam; pindelbam; checkbams; shimmerbam;}

process msi {
  executor 'local'
  publishDir "$params.output/$pid/dna_$params.projectid", mode: 'copy'
  errorStrategy 'ignore'
  input:
  set pid,tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from msibam
  output:
  file("${pid}*") into msiout
  script:
  """
  bash $baseDir/process_scripts/variants/msisensor.sh -r ${index_path} -p $pid -b $tumor -n $normal -c $capturebed
  """
}

process checkmates {
  executor 'local'
  publishDir "$params.output/$pid/dna_$params.projectid", mode: 'copy'
  errorStrategy 'ignore'
  input:
  set pid,tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from checkbams
  file(conf) from ncmconf
  output:
  file("${pid}*") into checkmateout
  script:
  """
  source /etc/profile.d/modules.sh
  module load python/2.7.x-anaconda git/v2.5.3
  python /project/shared/bicf_workflow_ref/seqprg/NGSCheckMate/ncm.py -B -d ./ -bed ${index_path}/NGSCheckMate.bed -O ./ -N ${pid}
  perl $baseDir/scripts/sequenceqc_somatic.pl -r ${index_path} -i ${pid}_all.txt -o ${pid}_${params.projectid}.sequence.stats.txt
  """
}
process sv {
  queue '32GB'
  errorStrategy 'ignore'
  publishDir "$params.output/$pid/dna_$params.projectid", mode: 'copy'

  input:
  set pid,tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from dellybam
  each algo from svalgo
  output:
  set pid,file("${pid}.${algo}.vcf.gz") into svvcf
  file("${pid}.${algo}.genefusion.txt") into svgf
  script:				       
  """
  bash $baseDir/process_scripts/variants/svcalling.sh -r $index_path -x ${tid} -y ${nid} -b $tumor -n $normal -p $pid -a ${algo} -f 
  """
}
process pindel {
  queue '128GB,256GB,256GBv1'
  errorStrategy 'ignore'
  publishDir "$params.output/$pid/dna_$params.projectid", mode: 'copy'
  input:
  set pid,tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from pindelbam
  output:
  file("${pid}.pindel_tandemdup.vcf.gz") into tdvcf
  set pid,file("${pid}.pindel.vcf.gz") into pindelvcf
  file("${pid}.pindel.genefusion.txt") into pindelgf
  script:
  """
  source /etc/profile.d/modules.sh
  module load samtools/gcc/1.8 snpeff/4.3q htslib/gcc/1.8 
  bash $baseDir/process_scripts/variants/svcalling.sh -r $index_path -p $pid -l ${index_path}/itd_genes.bed -a pindel -f
  """
}
process fb {
  queue '32GB'
  errorStrategy 'ignore'
  publishDir "$params.output/$pid/dna_$params.projectid", mode: 'copy'
  
  input:
  set pid,tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from fbbam
  each algo from fpsalgo
  output:
  set pid,file("${pid}.${algo}.vcf.gz") into fbvcf
  set pid,file("${pid}.${algo}.ori.vcf.gz") into fbori
  script:
  """
  bash $baseDir/process_scripts/variants/germline_vc.sh -r $index_path -p $pid -a ${algo}
  bash $baseDir/process_scripts/variants/uni_norm_annot.sh -g $snpeff_vers -r $index_path -p ${pid}.${algo} -v ${pid}.${algo}.vcf.gz
  """
}

process mutect {
  queue '128GB,256GB,256GBv1'
  errorStrategy 'ignore'
  publishDir "$params.output/$pid/dna_$params.projectid", mode: 'copy'
  input:
  set pid,tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from mutectbam
  output:
  set pid,file("${pid}.mutect.vcf.gz") into mutectvcf
  set pid,file("${pid}.mutect.ori.vcf.gz") into mutectori
  script:
  """
  bash $baseDir/process_scripts/variants/somatic_vc.sh $ponopt -r $index_path -p $pid -x $tid -y $nid -n $normal -t $tumor -a mutect
  bash $baseDir/process_scripts/variants/uni_norm_annot.sh -g $snpeff_vers -r $index_path -p ${pid}.mutect -v ${pid}.mutect.vcf.gz
  """
}
process ss {
  queue '32GB'
  errorStrategy 'ignore'
  publishDir "$params.output/$pid/dna_$params.projectid", mode: 'copy'
  input:
  set pid,tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from strelkabam
  each algo from ssalgo
  output:
  set pid,file("${pid}.${algo}.vcf.gz") into ssvcf
  set pid,file("${pid}.${algo}.ori.vcf.gz") into ssori
  script:
  """
  bash $baseDir/process_scripts/variants/somatic_vc.sh -r $index_path -p $pid -x $tid -y $nid -n $normal -t $tumor -a ${algo}
  bash $baseDir/process_scripts/variants/uni_norm_annot.sh -g $snpeff_vers -r $index_path -p ${pid}.${algo} -v ${pid}.${algo}.vcf.gz
  """
}

Channel
  .empty()
  .mix(mutectvcf,fbvcf,ssvcf,pindelvcf)
  .groupTuple(by:0)
  .set { vcflist}

process integrate {
  executor 'local'
  errorStrategy 'ignore'
  publishDir "$params.output/$pid", mode: 'copy'
  input:
  set pid,file(vcf) from vcflist
  file 'design.txt' from design_file
  output:
  file("${pid}_${params.projectid}.dna.vcf.gz") into unionvcf
  script:
  """
  source /etc/profile.d/modules.sh
  module load htslib/gcc/1.8
  bash $baseDir/process_scripts/variants/union.sh -r $index_path -p $pid
  cp ${pid}.union.vcf.gz ${pid}_${params.projectid}.dna.vcf.gz
  """
}
