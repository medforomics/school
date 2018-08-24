#!/usr/bin/env nextflow

params.input = './fastq'
params.output = './analysis'

params.fastqs="$params.input/*.fastq.gz"
params.design="$params.input/design.txt"

params.genome="/project/shared/bicf_workflow_ref/GRCh38"
params.capture="$params.genome/clinseq_prj/UTSWV2.bed"
params.pairs="pe"
params.cancer="detect"
params.markdups='picard'

reffa=file("$params.genome/genome.fa")
dbsnp="$params.genome/dbSnp.vcf.gz"
indel="$params.genome/GoldIndels.vcf.gz"

fastqs=file(params.fastqs)
design_file = file(params.design)
dbsnp=file(dbsnp)
knownindel=file(indel)
index_path = file(params.genome)
capture_bed = file(params.capture)

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
    iidx = header.findIndexOf{it == 'SampleID'};
    nidx = header.findIndexOf{it == 'SampleName'};
    fidx = header.findIndexOf{it == 'FamilyID'};
    oneidx = header.findIndexOf{it == 'FqR1'};
    twoidx = header.findIndexOf{it == 'FqR2'};
    if (twoidx == -1) {
       twoidx = oneidx
       }
    while (line = reader.readLine()) {
    	   def row = line.split("\t")
    if (fileMap.get(row[oneidx]) != null) {
	read << tuple(row[fidx],row[iidx],row[nidx],fileMap.get(row[oneidx]),fileMap.get(row[twoidx]))
	   }

}
}
if( ! read) { error "Didn't match any input files with entries in the design file" }

process trim {
  errorStrategy 'ignore'
  input:
  set subjid,sname,pair_id, file(read1), file(read2) from read
  output:
  set subjid,sname,pair_id, file("${pair_id}.trim.R1.fastq.gz"),file("${pair_id}.trim.R2.fastq.gz") into trimread
  file("${pair_id}.trimreport.txt") into trimstat 
  script:
  """
  bash $baseDir/process_scripts/preproc_fastq/trimgalore.sh -p ${pair_id} -a ${read1} -b ${read2}
  perl $baseDir/process_scripts/preproc_fastq/parse_trimreport.pl ${pair_id}.trimreport.txt *trimming_report.txt
  """
}
process align {
  errorStrategy 'ignore'
  //publishDir "$params.output", mode: 'copy'

  input:
  set subjid,sname,pair_id, file(fq1), file(fq2) from trimread
  output:
  set subjid,sname, file("${pair_id}.dedup.bam") into aligned
  file("${pair_id}.libcomplex.txt") into libcomplex
  """
  source /etc/profile.d/modules.sh
  touch ${pair_id}.dedup.stat.txt
  bash $baseDir/process_scripts/alignment/dnaseqalign.sh -r $index_path -p $pair_id -x $fq1 -y $fq2
  bash $baseDir/process_scripts/alignment/markdups.sh -a $params.markdups -b ${pair_id}.bam -p $pair_id
  mv ${pair_id}.dedup.stat.txt ${pair_id}.libcomplex.txt
  """
}

aligned
   .groupTuple(by:0)
   .set {bamgrp}

process mergebam {
  publishDir "$params.output/$subjid/$pair_id", mode: 'copy'

  input:
  set subjid,pair_id,file(bams) from bamgrp
  output:
  set subjid,pair_id,file("${pair_id}.bam"),file("${pair_id}.bam.bai") into qcbam
  set subjid,pair_id,file("${pair_id}.bam"),file("${pair_id}.bam.bai") into qcbam2
  set subjid,pair_id,file("${pair_id}.bam"),file("${pair_id}.bam.bai") into deduped

  script:
  """
  source /etc/profile.d/modules.sh
  module load samtools/1.6
  count=\$(ls *.bam |wc -l)
  if [ \$count -gt 1 ]
  then
  samtools merge -@ \$SLURM_CPUS_ON_NODE merge.bam *.bam
  else
  mv *.bam merge.bam
  fi
  samtools sort --threads \$SLURM_CPUS_ON_NODE -o ${pair_id}.bam merge.bam
  samtools index ${pair_id}.bam
  bash $baseDir/process_scripts/alignment/bam2tdf.sh -r $index_path -b ${pair_id}.bam -p ${pair_id}.raw
  """
}
process uniqqc {
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid/$pair_id", mode: 'copy'

  input:
  set subjid,pair_id, file(sbam),file(sbai) from qcbam2
  output:
  file("${pair_id}.dedupcov.txt") into dedupcov
  file("${pair_id}.covuniqhist.txt") into covuniqhist
  file("*coverageuniq.txt") into covuniqstat
  file("${pair_id}.uniq.tdf") into uniqtdf
  script:
  """
  module load samtools/1.6
  samtools view -1 -F 1024 -o ${pair_id}.rmdup.bam $sbam
  bash $baseDir/process_scripts/alignment/bam2tdf.sh -r $index_path -b ${pair_id}.rmdup.bam -p ${pair_id}.uniq
  bash $baseDir/process_scripts/alignment/bamqc.sh -c $capture_bed -n dna -r $index_path -b ${pair_id}.rmdup.bam -p $pair_id
  mv ${pair_id}.genomecov.txt ${pair_id}.dedupcov.txt
  mv ${pair_id}.covhist.txt ${pair_id}.covuniqhist.txt
  mv ${pair_id}_lowcoverage.txt ${pair_id}_lowcoverageuniq.txt
  mv ${pair_id}_exoncoverage.txt ${pair_id}_exoncoverageuniq.txt
  """
}
process seqqc {
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid/$pair_id", mode: 'copy'

  input:
  set subjid,pair_id, file(sbam),file(sbai) from qcbam
  output:
  file("${pair_id}.flagstat.txt") into alignstats
  file("${pair_id}.ontarget.flagstat.txt") into ontarget
  file("${pair_id}.meanmap.txt") into meanmap
  file("${pair_id}.hist.txt") into insertsize
  file("${pair_id}.alignmentsummarymetrics.txt") into alignmentsummarymetrics
  file("*fastqc*") into fastqc
  set pair_id, file("${pair_id}.ontarget.bam"),file("${pair_id}.ontarget.bam.bai") into ontargetbam
  set pair_id,file("${pair_id}.ontarget.bam"),file("${pair_id}.ontarget.bam.bai") into genocovbam
  file("${pair_id}.genomecov.txt") into genomecov
  file("${pair_id}.covhist.txt") into covhist
  file("*coverage.txt") into capcovstat
  file("${pair_id}.mapqualcov.txt") into mapqualcov
  script:
  """
  bash $baseDir/process_scripts/alignment/bamqc.sh -c $capture_bed -n dna -r $index_path -b $sbam -p $pair_id
  """
}
process parse_stat {
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'

  input:
  file(txt) from alignstats.toList()
  file(lc) from libcomplex.toList()
  file(is) from insertsize.toList()
  file(gc) from genomecov.toList()
  file(on) from ontarget.toList()
  file(tr) from trimstat.toList()
  file(mq) from mapqualcov.toList()
  file(de) from dedupcov.toList()
  file(mm) from meanmap.toList()
  file(asmet) from alignmentsummarymetrics.toList()
  
  output:
  file('*sequence.stats.txt')
  file('*.png')
  script:
  """
  source /etc/profile.d/modules.sh
  module load R/3.2.1-intel git/gcc/v2.12.2
  perl $baseDir/scripts/sequenceqc_alignment.pl -r ${index_path} *.genomecov.txt
  perl $baseDir/scripts/covstats.pl *.mapqualcov.txt
  Rscript $baseDir/scripts/plot_hist_genocov.R
  """
}

process gatkbam {
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid/$pair_id", mode: 'copy'

  input:
  set subjid,pair_id, file(sbam),file(sbai) from deduped

  output:
  set subjid,pair_id,file("${pair_id}.final.bam"),file("${pair_id}.final.bai") into gatkbam
  
  script:
  """
  bash $baseDir/process_scripts/variants/gatkrunner.sh -a gatkbam -b $sbam -r ${index_path} -p $pair_id
  """
}
