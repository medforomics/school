#!/usr/bin/env nextflow

// Default parameter values to run tests
params.input = './'
params.output = './'

params.bams="$params.input/*.bam"
params.design="$params.input/design.txt"
params.callsvs="skip"
params.genome="/project/shared/bicf_workflow_ref/GRCh38"
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
    totidx = header.findIndexOf{it == 'TumorFinalBAM'};
    notidx = header.findIndexOf{it == 'NormalFinalBAM'};

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
  errorStrategy 'ignore'
  input:
  set pid,tid,nid,file(tumor),file(normal) from oribam
  output:
  set pid,tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into dellybam
  set pid,tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into mantrabam
  set pid,tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into checkbams

  script:
  """
  bash $baseDir/process_scripts/alignment/indexbams.sh 
  """
}
process indextarbams {
  errorStrategy 'ignore'
  input:
  set pid,tid,nid,file(tumor),file(normal) from tarbam
  output:
  set pid,tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into mutectbam
  set pid,tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into strelkabam
  set pid,tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into ssbam
  set pid,tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into shimmerbam
  set pid,tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into vscanbam
  set pid,tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into virmidbam
  set pid,tid,nid into pairnames
  script:
  """
  bash $baseDir/process_scripts/alignment/indexbams.sh 
  """
}
process checkmates {
  publishDir "$params.output/$pid/somatic", mode: 'copy'
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

process delly {
  publishDir "$params.output/$pid/somatic", mode: 'copy'
  errorStrategy 'ignore'
  input:
  set pid,tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from dellybam

  output:
  file("${pid}.delly.vcf.gz") into dellyvcf
  when:
  params.callsvs == "detect"
  script:
  """
  source /etc/profile.d/modules.sh
  perl $baseDir/scripts/make_delly_sample.pl ${tid} ${nid}
  bash $baseDir/process_scripts/variants/svcalling.sh -r ${index_path} -p ${pid} -b ${tumor} -n ${normal} -i ${tid} -m delly
  """
}

process sstumor {
  errorStrategy 'ignore'
  publishDir "$params.output/$pid/somatic", mode: 'copy'
  input:
  set pid,tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from ssbam
    
  output:
  set pid, file("${pid}.sssom.vcf.gz") into ssvcf
  set pid,file("${pid}.sssom.annot.vcf.gz") into ssannot
  set pid,file("${pid}.sssom.ori.vcf.gz") into ssori

  script:
  """
  bash $baseDir/process_scripts/variants/somatic_vc.sh -r $index_path -p $pid -x $tid -y $nid -n $normal -t $tumor -a speedseq
  bash $baseDir/process_scripts/variants/norm_annot.sh -r $index_path -p ${pid}.sssom -v ${pid}.sssom.vcf.gz
  mv ${pid}.sssom.vcf.gz ${pid}.sssom.ori.vcf.gz
  mv ${pid}.sssom.norm.vcf.gz ${pid}.sssom.vcf.gz
  bash $baseDir/process_scripts/variants/annotvcf.sh -p ${pid}.sssom -r $index_path -v ${pid}.sssom.vcf.gz

  """
}
process mutect {
  errorStrategy 'ignore'
  publishDir "$params.output/$pid/somatic", mode: 'copy'

  input:
  set pid,tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from mutectbam

  output:
  set pid,file("${pid}.mutect.vcf.gz") into mutectvcf
  set pid,file("${pid}.mutect.ori.vcf.gz") into mutectori
  set pid,file("${pid}.mutect.annot.vcf.gz") into mutectannot
  script:
  """
  bash $baseDir/process_scripts/variants/somatic_vc.sh -r $index_path -p $pid -x $tid -y $nid -n $normal -t $tumor -a mutect2
  bash $baseDir/process_scripts/variants/norm_annot.sh -r $index_path -p ${pid}.mutect -v ${pid}.mutect.vcf.gz
  mv ${pid}.mutect.vcf.gz ${pid}.mutect.ori.vcf.gz
  mv ${pid}.mutect.norm.vcf.gz ${pid}.mutect.vcf.gz
  bash $baseDir/process_scripts/variants/annotvcf.sh -p ${pid}.mutect -r $index_path -v ${pid}.mutect.vcf.gz
  """
}
// Channel
//   .empty()
//   .mix(mantrabam,strelkabam)
//   .groupTuple(by:0)
//   .into { illuminabams }

// process strelka {
//   errorStrategy 'ignore'
//   publishDir "$params.output/$pid/somatic", mode: 'copy'

//   input:
//   set pid,mtid,mnid,file(mtumor),file(mnormal),file(mtidx),file(mnidx),tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from illuminabams

//   output:
//   set pid,file("${pid}.strelka2.vcf.gz") into strelkavcf
//   script:
//   """
//   bash $baseDir/process_scripts/variants/somatic_vc.sh -r $index_path -p $pid -x $tid -y $nid -n $normal -t $tumor -a strelka2
//   """
// }
process varscan {
  errorStrategy 'ignore'
  publishDir "$params.output/$pid/somatic", mode: 'copy'
  input:
  set pid,tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from vscanbam
  output:
  set pid,file("${pid}.varscan.vcf.gz") into varscanvcf
  set pid,file("${pid}.varscan.ori.vcf.gz") into varscanori
  set pid,file("${pid}.varscan.annot.vcf.gz") into varscannot
  set pid,file("${pid}.vscancnv.copynumber.txt") into varscancnv
  script:
  """
  source /etc/profile.d/modules.sh
  bash $baseDir/process_scripts/variants/somatic_vc.sh -r $index_path -p $pid -x $tid -y $nid -n $normal -t $tumor -a varscan
  mv vscancnv.copynumber ${pid}.vscancnv.copynumber.txt
  bash $baseDir/process_scripts/variants/norm_annot.sh -r $index_path -p ${pid}.varscan -v ${pid}.varscan.vcf.gz
  mv ${pid}.varscan.vcf.gz ${pid}.varscan.ori.vcf.gz
  mv ${pid}.varscan.norm.vcf.gz ${pid}.varscan.vcf.gz
  bash $baseDir/process_scripts/variants/annotvcf.sh -p ${pid}.varscan -r $index_path -v ${pid}.varscan.vcf.gz
  """
}

process shimmer {
  errorStrategy 'ignore'
  publishDir "$params.output/$pid/somatic", mode: 'copy'
  input:
  set pid,tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from shimmerbam
  output:
  set pid, file("${pid}.shimmer.vcf.gz") into shimmervcf
  set pid, file("${pid}.shimmer.vcf.gz") into shimmerori
  set pid, file("${pid}.shimmer.vcf.gz") into shimmerannot
  script:
  """
  bash $baseDir/process_scripts/variants/somatic_vc.sh -r $index_path -p $pid -x $tid -y $nid -n $normal -t $tumor -a shimmer
  bash $baseDir/process_scripts/variants/norm_annot.sh -r $index_path -p ${pid}.shimmer -v ${pid}.shimmer.vcf.gz
  mv ${pid}.shimmer.vcf.gz ${pid}.shimmer.ori.vcf.gz
  mv ${pid}.shimmer.norm.vcf.gz ${pid}.shimmer.vcf.gz
  bash $baseDir/process_scripts/variants/annotvcf.sh -p ${pid}.shimmer -r $index_path -v ${pid}.shimmer.vcf.gz
  """
}

process virmid {
  errorStrategy 'ignore'
  publishDir "$params.output/$pid/somatic", mode: 'copy'
  input:
  set pid,tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from virmidbam
  output:
  set pid, file("${pid}.virmid.vcf.gz") into virmidvcf
  set pid, file("${pid}.virmid.annot.vcf.gz") into virmidannot
  set pid, file("${pid}.virmid.ori.vcf.gz") into virmidori
  script:
  """
  bash $baseDir/process_scripts/variants/somatic_vc.sh -r $index_path -p $pid -x $tid -y $nid -n $normal -t $tumor -a virmid
  bash $baseDir/process_scripts/variants/norm_annot.sh -r $index_path -p ${pid}.virmid -v ${pid}.virmid.vcf.gz
  mv ${pid}.virmid.vcf.gz ${pid}.virmid.ori.vcf.gz
  mv ${pid}.virmid.norm.vcf.gz ${pid}.virmid.vcf.gz
  bash $baseDir/process_scripts/variants/annotvcf.sh -p ${pid}.virmid -r $index_path -v ${pid}.virmid.vcf.gz
  """
}

Channel
  .empty()
  .mix(ssvcf,mutectvcf,varscanvcf,virmidvcf,shimmervcf)
  .groupTuple(by:0)
  .set { vcflist}


process integrate {
  errorStrategy 'ignore'
  publishDir "$params.output/$pid/somatic", mode: 'copy'
  input:
  set pid,file(vcf) from vcflist
  file 'design.txt' from design_file
  output:
  file("${pid}${params.projectid}.somaticunion.vcf.gz") into union
  file("${pid}${params.projectid}.somatic.vcf.gz") into annotvcf
  script:
  """
  source /etc/profile.d/modules.sh
  module load samtools/1.6
  bash $baseDir/process_scripts/variants/union.sh -r $index_path -p $pid
  bash $baseDir/process_scripts/variants/annotvcf.sh -p $pid -r $index_path -v ${pid}.union.vcf.gz
  perl $baseDir/scripts/somatic_filter.pl ${pid}.annot.vcf.gz
  bgzip ${pid}.somatic.vcf
  mv ${pid}.somatic.vcf.gz ${pid}${params.projectid}.somatic.vcf.gz
  mv ${pid}.union.vcf.gz ${pid}${params.projectid}.somaticunion.vcf.gz
  """
}
