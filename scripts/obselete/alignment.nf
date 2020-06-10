#!/usr/bin/env nextflow

params.input = './fastq'
params.output = './analysis'

params.fastqs="$params.input/*.fastq.gz"
params.design="$params.input/design.txt"

params.genome="/project/shared/bicf_workflow_ref/human/grch38_cloud/dnaref"
params.virus_genome="/project/shared/bicf_workflow_ref/human_virus_genome/clinlab_idt_genomes"

params.pairs="pe"
params.markdups='fgbio_umi'
params.aligner='bwa'
params.cnv = "detect"

params.version = 'v4'

reffa=file("$params.genome/genome.fa")
dbsnp="$params.genome/dbSnp.vcf.gz"
indel="$params.genome/GoldIndels.vcf.gz"

fastqs=file(params.fastqs)
design_file=file(params.design)
dbsnp=file(dbsnp)
knownindel=file(indel)
index_path=file(params.genome)
capture_bed = file(params.capture)
capturedir = file(params.capturedir)
virus_index_path=file(params.virus_genome)

skipCNV = false
if(capturedir.isEmpty()) {
  skipCNV = true
}

alignopts = ''
if (params.markdups == 'fgbio_umi') {
   alignopts='-u'
}

def fileMap = [:]

fastqs.each {
    final fileName = it.getFileName().toString()
    prefix = fileName.lastIndexOf('/')
    fileMap[fileName] = it
}
def read = []
new File(params.design).withReader { reader ->
    def hline = reader.readLine()
    def header = hline.split("\t")
    prefixidx = header.findIndexOf{it == 'SampleID'};
    oneidx = header.findIndexOf{it == 'FqR1'};
    twoidx = header.findIndexOf{it == 'FqR2'};
    fidx = header.findIndexOf{it == 'FamilyID'};
    sidx = header.findIndexOf{it == 'SubjectID'};
    if (sidx == -1) {
       sidx = prefixidx
    }
    if (fidx == -1) {
       fidx = sidx
    }
    if (twoidx == -1) {
       twoidx = oneidx
       }
    while (line = reader.readLine()) {
    	   def row = line.split("\t")
    if (fileMap.get(row[oneidx]) != null) {
	read << tuple(row[fidx],row[prefixidx],fileMap.get(row[oneidx]),fileMap.get(row[twoidx]))
	   }

}
}
if( ! read) { error "Didn't match any input files with entries in the design file" }

process dtrim {
  queue '32GB,128GB,256GB,256GBv1'
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid/$pair_id", mode: 'copy'
  input:
  set subjid,pair_id, file(read1), file(read2) from read
  output:
  set subjid,pair_id, file("${pair_id}.trim.R1.fastq.gz"),file("${pair_id}.trim.R2.fastq.gz"),file("${pair_id}.trimreport.txt") into trimread
  script:
  """
  bash $baseDir/process_scripts/preproc_fastq/trimgalore.sh -p ${pair_id} -a ${read1} -b ${read2} -f
  """
}

process dalign {
  queue '32GB,128GB,256GB,256GBv1'
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid/$pair_id", mode: 'copy'

  input:
  set subjid,pair_id, file(fq1), file(fq2), file(trimreport) from trimread
  output:
  set subjid,pair_id, file("${pair_id}.bam") into aligned
  set subjid,pair_id, file("${pair_id}.bam"),file("${pair_id}.bam.bai"),file(trimreport) into aligned2
  set subjid,pair_id, file("${pair_id}.bam") into virusalign
  set subjid,pair_id, file("${pair_id}.bam"), file("${pair_id}.bam.bai") into cnvbam
  set subjid,pair_id, file("${pair_id}.bam"), file("${pair_id}.bam.bai") into itdbam

  """
  bash $baseDir/process_scripts/alignment/dnaseqalign.sh -r $index_path -p $pair_id -x $fq1 -y $fq2 $alignopts
  """
 }

process valign {
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid/$pair_id", mode: 'copy'

  input:
  set subjid,pair_id, file(sbam) from virusalign
  output:
  file("${pair_id}.viral.seqstats.txt") into viralseqstats
  """
  bash $baseDir/process_scripts/alignment/virusalign.sh -b ${pair_id}.bam -p ${pair_id} -r $virus_index_path -f
  """
}

process markdups_consensus {
  errorStrategy 'ignore'
  queue '32GB,128GB,256GB,256GBv1'
  publishDir "$params.output/$subjid/$pair_id", mode: 'copy'

  input:
  set subjid, pair_id, file(sbam) from aligned
  output:
  set subjid, pair_id, file("${pair_id}.consensus.bam"),file("${pair_id}.consensus.bam.bai") into togatkbam
  file("*.txt") into statfile
  script:
  """
  bash $baseDir/process_scripts/alignment/markdups.sh -a $params.markdups -b $sbam -p $pair_id -r $index_path
  mv ${pair_id}.dedup.bam ${pair_id}.consensus.bam
  mv ${pair_id}.dedup.bam.bai ${pair_id}.consensus.bam.bai
  if [[ ! -f ${pair_id}.group.bam ]] 
  then
    cp ${pair_id}.consensus.bam ${pair_id}.group.bam
    cp ${pair_id}.consensus.bam.bai ${pair_id}.group.bam.bai
  fi
  """
}

process qc_gbam           {
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid/$pair_id", mode: 'copy'
  queue '128GB,256GB,256GBv1'

  input:
  set subjid,pair_id, file(sbam),file(bai),file(trimreport) from aligned2
  output:
  file("*fastqc*") into fastqc
  file("${pair_id}*.txt") into alignstats
  script:
  """
  bash $baseDir/process_scripts/alignment/bamqc.sh -c $capture_bed -n dna -r $index_path -b ${sbam} -p $pair_id  
  perl $baseDir/scripts/sequenceqc_dna.pl -r ${index_path} -e $params.version *.genomecov.txt
  """
}

process cnv {
  queue '32GB,128GB,256GB,256GBv1'
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid/$pair_id", mode: 'copy'
  input:
  set subjid,pair_id,file(sbam),file(sidx) from cnvbam

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
  source /etc/profile.d/modules.sh
  module load htslib/gcc/1.8 
  bash $baseDir/process_scripts/variants/cnvkit.sh -r $index_path -b $sbam -p $pair_id -d $capturedir
  """
}

process itdseek {
  queue '32GB,128GB,256GB,256GBv1'
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid/$pair_id", mode: 'copy'
  input:
  set subjid,pair_id,file(sbam),file(sidx) from itdbam

  output:
  file("${pair_id}.itdseek_tandemdup.vcf.gz") into itdseekvcf

  script:
  """
  bash $baseDir/process_scripts/variants/svcalling.sh -b $sbam -r $index_path -p $pair_id -l ${index_path}/itd_genes.bed -a itdseek -f
  """
}

process gatkbam {
  queue '32GB,128GB,256GB,256GBv1'
  publishDir "$params.output/$subjid/$pair_id", mode: 'copy'

  input:
  set subjid, pair_id, file(sbam),file(idx) from togatkbam
  output:
  set subjid, pair_id, file("${pair_id}.final.bam*") into gatkbam
  script:
  """
  bash $baseDir/process_scripts/variants/gatkrunner.sh -a gatkbam -b $sbam -r $index_path -p $pair_id
  """
}
