#!/usr/bin/env nextflow

params.input = './fastq'
params.output = './analysis'

params.bams="$params.input/*.bam"
params.design="$params.input/design.txt"

params.genome="/project/shared/bicf_workflow_ref/GRCh38"
params.capture="$params.genome/clinseq_prj/UTSWV2.bed"
params.pairs="pe"
params.markdups='fgbio_umi'

reffa=file("$params.genome/genome.fa")
dbsnp="$params.genome/dbSnp.vcf.gz"
indel="$params.genome/GoldIndels.vcf.gz"

bams=file(params.bams)
design_file = file(params.design)
dbsnp=file(dbsnp)
knownindel=file(indel)
index_path = file(params.genome)
capture_bed = file(params.capture)

alignopts = ''
if (params.markdups == 'fgbio_umi') {
   alignopts='-u'
}

def fileMap = [:]

bams.each {
    final fileName = it.getFileName().toString()
    prefix = fileName.lastIndexOf('/')
    fileMap[fileName] = it
}
def oribam = []
new File(params.design).withReader { reader ->
    def hline = reader.readLine()
    def header = hline.split("\t")
    tidx = header.findIndexOf{it == 'SampleID'};
    oneidx = header.findIndexOf{it == 'BAM'};
    fidx = header.findIndexOf{it == 'FamilyID'};
    sidx = header.findIndexOf{it == 'SubjectID'};
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
	   }

}
}
process indexoribams {
  errorStrategy 'ignore'
  input:
  set sid,tid,file(tumor) from oribam
  output:
  set sid,tid,file(tumor),file("${tumor}.bai") into aligned
  set sid,tid,file(tumor),file("${tumor}.bai") into aligned2
  script:
  """
  bash $baseDir/process_scripts/alignment/indexbams.sh 
  """
}

process markdups_consensus {
  publishDir "$params.output/$subjid/$pair_id", mode: 'copy'

  input:
  set subjid, pair_id, file(sbam) from aligned
  output:
  set subjid, pair_id, file("${pair_id}.consensus.bam") into deduped
  set subjid, pair_id, file("${pair_id}.consensus.bam") into qcbam
  script:
  """
  bash $baseDir/process_scripts/alignment/markdups.sh -a fgbio_umi -b $sbam -p $pair_id
  mv ${pair_id}.dedup.bam ${pair_id}.consensus.bam
  """
}
process markdups_picard {
  publishDir "$params.output/$subjid/$pair_id", mode: 'copy'

  input:
  set subjid, pair_id, file(sbam) from aligned2
  output:
  file("*fastqc*") into fastqc
  file("${pair_id}.flagstat.txt") into alignstats
  file("${pair_id}.ontarget.flagstat.txt") into ontarget
  file("${pair_id}.meanmap.txt") into meanmap
  file("${pair_id}.libcomplex.txt") into libcomplex
  file("${pair_id}.hist.txt") into insertsize
  file("${pair_id}.alignmentsummarymetrics.txt") into alignmentsummarymetrics
  file("${pair_id}.genomecov.txt") into genomecov
  file("${pair_id}.covhist.txt") into covhist
  file("*coverage.txt") into capcovstat
  file("${pair_id}.mapqualcov.txt") into mapqualcov

  script:
  """
  bash $baseDir/process_scripts/alignment/markdups.sh -a picard_umi -b $sbam -p $pair_id
  bash $baseDir/process_scripts/alignment/bamqc.sh -c $capture_bed -n dna -r $index_path -b ${pair_id}.dedup.bam -p $pair_id  
  """
}

process seqqc {
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid/$pair_id", mode: 'copy'

  input:
  set subjid, pair_id, file(sbam) from qcbam
  output:
  file("${pair_id}.dedupcov.txt") into dedupcov
  file("${pair_id}.covuniqhist.txt") into covuniqhist
  file("*coverageuniq.txt") into covuniqstat
  script:
  """
  bash $baseDir/process_scripts/alignment/bamqc.sh -c $capture_bed -n dna -r $index_path -b $sbam -p $pair_id
  mv ${pair_id}.genomecov.txt ${pair_id}.dedupcov.txt
  mv ${pair_id}.covhist.txt ${pair_id}.covuniqhist.txt
  mv ${pair_id}_lowcoverage.txt ${pair_id}_lowcoverageuniq.txt
  mv ${pair_id}_exoncoverage.txt ${pair_id}_exoncoverageuniq.txt
  """
}

process gatkbam {
  //errorStrategy 'ignore'
  publishDir "$params.output/$subjid/$pair_id", mode: 'copy'

  input:
  set subjid, pair_id, file(sbam) from deduped

  output:
  set subjid, pair_id,file("${pair_id}.final.bam"),file("${pair_id}.final.bai") into gatkbam
  
  script:
  """
  bash $baseDir/process_scripts/variants/gatkrunner.sh -a gatkbam -b $sbam -r ${index_path} -p $pair_id
  """
}
