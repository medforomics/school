#!/usr/bin/env nextflow

// Default parameter values to run tests
params.input = './'
params.output = './'

params.bams="$params.input/*.bam"
params.design="$params.input/design.txt"
params.callsvs="skip"
params.genome="/project/shared/bicf_workflow_ref/human/GRCh38"
params.projectid=''

dbsnp="$params.genome/dbSnp.vcf.gz"
cosmic="$params.genome/cosmic.vcf.gz"

design_file = file(params.design)
bams=file(params.bams)

index_path = file(params.genome)
ncmconf = file("$params.genome/ncm.conf")
dbsnp=file(dbsnp)

snpeff_vers = 'GRCh38.86';

def fileMap = [:]

bams.each {
    final fileName = it.getFileName().toString()
    prefix = fileName.lastIndexOf('/')
    fileMap[fileName] = it
}
def oribam = []
def tarbam = []
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
  set pid,tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into pindelbam
  set pid,tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into platbam
  set pid,tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into fbbam
  set pid,tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into checkbams
  set pid,tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into strelkabam
  set pid,tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into shimmerbam
  set pid,tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into virmidbam
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
process checkmates {
  executor 'local'
  publishDir "$params.output/$pid/somatic$params.projectid", mode: 'copy'
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
  perl $baseDir/scripts/sequenceqc_somatic.pl -r ${index_path} -i ${pid}_all.txt -o ${pid}${params.projectid}.sequence.stats.txt
  """
}
process pindel {
  queue '32GB,128GB,256GB,256GBv1'
  errorStrategy 'ignore'
  publishDir "$params.output/$pid", mode: 'copy'
  input:
  set pid,tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from pindelbam
  output:
  file("${pid}.pindel_tandemdup.pass.vcf.gz") into tdvcf
  file("${pid}.pindel_indel.pass.vcf.gz") into pindelvcf
  file("${pid}.dna.genefusion.txt") into gf
  script:
  """
  source /etc/profile.d/modules.sh
  bash $baseDir/process_scripts/variants/pindel.sh -r ${index_path} -p ${pid}
  perl $baseDir/process_scripts/variants/filter_pindel.pl -d ${pid}.pindel_tandemdup.vcf.gz -s ${pid}.pindel_sv.vcf.gz -i ${pid}.pindel_indel.vcf.gz
  module load htslib/gcc/1.8 snpeff/4.3q
  bgzip ${pid}.pindel_indel.pass.vcf
  bgzip ${pid}.pindel_tandemdup.pass.vcf
  grep '#CHROM' ${pid}.pindel_sv.pass.vcf > ${pid}.dna.genefusion.txt
  cat ${pid}.pindel_sv.pass.vcf | \$SNPEFF_HOME/scripts/vcfEffOnePerLine.pl |java -jar \$SNPEFF_HOME/SnpSift.jar extractFields - CHROM POS END ANN[*].EFFECT ANN[*].GENE ANN[*].HGVS_C ANN[*].HGVS_P GEN[*] |grep -E 'CHROM|gene_fusion' |uniq >> ${pid}.dna.genefusion.txt
  """
}
process freebayes {
  queue '32GB'
  errorStrategy 'ignore'
  publishDir "$params.output/$pid/somatic$params.projectid", mode: 'copy'

  input:
  set pid,tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from fbbam
  output:
  set pid,file("${pid}.fb.vcf.gz") into fbvcf
  set pid,file("${pid}.fb.ori.vcf.gz") into fbori
  script:
  """
  bash $baseDir/process_scripts/variants/germline_vc.sh -r $index_path -p $pid -a freebayes
  bash $baseDir/process_scripts/variants/uni_norm_annot.sh -r $index_path -p ${pid}.fb -v ${pid}.freebayes.vcf.gz
  """
}
process platypus {
  queue '32GB'
  errorStrategy 'ignore'
  publishDir "$params.output/$pid/somatic$params.projectid", mode: 'copy'

  input:
  set pid,tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from platbam
  output:
  set pid,file("${pid}.platypus.vcf.gz") into platvcf
  set pid,file("${pid}.platypus.ori.vcf.gz") into platori
  script:				       
  """
  bash $baseDir/process_scripts/variants/germline_vc.sh -r $index_path -p $pid -a platypus
  bash $baseDir/process_scripts/variants/uni_norm_annot.sh -r $index_path -p ${pid}.platypus -v ${pid}.platypus.vcf.gz
  """
}

process mutect {
  queue '128GB,256GB,256GBv1'
  errorStrategy 'ignore'
  publishDir "$params.output/$pid/somatic$params.projectid", mode: 'copy'
  input:
  set pid,tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from mutectbam
  output:
  set pid,file("${pid}.mutect.vcf.gz") into mutectvcf
  set pid,file("${pid}.mutect.ori.vcf.gz") into mutectori
  script:
  """
  bash $baseDir/process_scripts/variants/somatic_vc.sh -r $index_path -p $pid -x $tid -y $nid -n $normal -t $tumor -a mutect2
  bash $baseDir/process_scripts/variants/uni_norm_annot.sh -r $index_path -p ${pid}.mutect -v ${pid}.mutect.vcf.gz
  """
}
process strelka {
  queue '32GB'
  errorStrategy 'ignore'
  publishDir "$params.output/$pid/somatic$params.projectid", mode: 'copy'
  input:
  set pid,tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from strelkabam
  output:
  set pid,file("${pid}.strelka2.vcf.gz") into strelkavcf
  set pid,file("${pid}.strelka2.ori.vcf.gz") into strelkaori
  script:
  """
  bash $baseDir/process_scripts/variants/somatic_vc.sh -r $index_path -p $pid -x $tid -y $nid -n $normal -t $tumor -a strelka2
  bash $baseDir/process_scripts/variants/uni_norm_annot.sh -r $index_path -p ${pid}.strelka2 -v ${pid}.strelka2.vcf.gz
  """
}
process shimmer {
  queue '32GB'
  errorStrategy 'ignore'
  publishDir "$params.output/$pid/somatic$params.projectid", mode: 'copy'
  input:
  set pid,tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from shimmerbam
  output:
  set pid, file("${pid}.shimmer.vcf.gz") into shimmervcf
  set pid, file("${pid}.shimmer.ori.vcf.gz") into shimmerori
  script:
  """
  bash $baseDir/process_scripts/variants/somatic_vc.sh -r $index_path -p $pid -x $tid -y $nid -n $normal -t $tumor -a shimmer
  bash $baseDir/process_scripts/variants/uni_norm_annot.sh -r $index_path -p ${pid}.shimmer -v ${pid}.shimmer.vcf.gz
  """
}

process virmid {
  queue '32GB'
  errorStrategy 'ignore'
  publishDir "$params.output/$pid/somatic$params.projectid", mode: 'copy'
  input:
  set pid,tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from virmidbam
  output:
  set pid, file("${pid}.virmid.vcf.gz") into virmidvcf
  set pid, file("${pid}.virmid.ori.vcf.gz") into virmidori
  script:
  """
  bash $baseDir/process_scripts/variants/somatic_vc.sh -r $index_path -p $pid -x $tid -y $nid -n $normal -t $tumor -a virmid
  bash $baseDir/process_scripts/variants/uni_norm_annot.sh -r $index_path -p ${pid}.virmid -v ${pid}.virmid.vcf.gz
  """
}

Channel
  .empty()
  .mix(mutectvcf,platvcf,fbvcf,shimmervcf,strelkavcf)
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
//  file("${pid}${params.projectid}.somatic.vcf.gz") into annotvcf
  file("${pid}${params.projectid}.dna.vcf.gz") into unionvcf
  script:
  """
  source /etc/profile.d/modules.sh
  module load htslib/gcc/1.8
  bash $baseDir/process_scripts/variants/union.sh -r $index_path -p $pid
  cp ${pid}.union.vcf.gz ${pid}_${params.projectid}.dna.vcf.gz
  #ln -s  ${pid}.union.vcf.gz ${pid}.annot.vcf.gz
  #perl $baseDir/scripts/somatic_filter.pl ${pid}.annot.vcf.gz
  #bgzip ${pid}.somatic.vcf
  #mv ${pid}.somatic.vcf.gz ${pid}${params.projectid}.somatic.vcf.gz
  """
}
