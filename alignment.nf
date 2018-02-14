#!/usr/bin/env nextflow

params.input = './fastq'
params.output = './analysis'

params.fastqs="$params.input/*.fastq.gz"
params.design="$params.input/design.txt"

params.genome="/project/shared/bicf_workflow_ref/GRCh38"
params.capture="$params.genome/UTSWV2.bed"
params.pairs="pe"
params.markdups='fgbio_umi'

reffa=file("$params.genome/genome.fa")
dbsnp="$params.genome/dbSnp.vcf.gz"
indel="$params.genome/GoldIndels.vcf.gz"

fastqs=file(params.fastqs)
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

fastqs.each {
    final fileName = it.getFileName().toString()
    prefix = fileName.lastIndexOf('/')
    fileMap[fileName] = it
}
def prefix = []
def read1_files = []
def read2_files = []
def mername = []
new File(params.design).withReader { reader ->
    def hline = reader.readLine()
    def header = hline.split("\t")
    prefixidx = header.findIndexOf{it == 'SampleID'};
    oneidx = header.findIndexOf{it == 'FqR1'};
    twoidx = header.findIndexOf{it == 'FqR2'};
    if (twoidx == -1) {
       twoidx = oneidx
       }
    while (line = reader.readLine()) {
    	   def row = line.split("\t")
    if (fileMap.get(row[oneidx]) != null) {
	prefix << row[prefixidx]
	read1_files << fileMap.get(row[oneidx])
        read2_files << fileMap.get(row[twoidx])
	   }

}
}
if( ! prefix) { error "Didn't match any input files with entries in the design file" }

Channel
  .from(read1_files)
  .set { read1 }
  
Channel
  .from(read2_files)
  .set { read2 }

Channel
  .from(prefix)
  .set { read_pe }

Channel
  .from(1..10)
  .set {counter}

process trimpe {
  errorStrategy 'ignore'
  //publishDir "$params.output/$pair_id", mode: 'copy'
  input:
  val pair_ids from read_pe.buffer(size: 20, remainder: true)
  file(read1s) from read1.buffer(size: 20, remainder: true)
  file(read2s) from read2.buffer(size: 20, remainder: true)
  val ct from counter
  output:
  file("*_val_1.fq.gz") into trimpe_r1s mode flatten
  file("*_val_2.fq.gz") into trimpe_r2s mode flatten
  file("*.trimreport.txt") into trimstat mode flatten

  script:
  assert pair_ids.size() == read1s.size()
  assert pair_ids.size() == read2s.size()
  def cmd = ''
  for( int i=0; i<pair_ids.size(); i++){
    cmd +="mv ${read1s[i]} ${pair_ids[i]}.R1.fq.gz\n"
    cmd +="mv ${read2s[i]} ${pair_ids[i]}.R2.fq.gz\n"
    cmd +="trim_galore --paired --stringency 3 -q 25 --illumina --gzip --length 35 ${pair_ids[i]}.R1.fq.gz ${pair_ids[i]}.R2.fq.gz &\n"
  }

  """
  source /etc/profile.d/modules.sh
  module load trimgalore/0.4.1 cutadapt/1.9.1
  ${cmd}
  wait
  perl $baseDir/scripts/parse_trimreport.pl trimreport_${ct}.txt *trimming_report.txt
  """
}
trimpe_r1s
  .map { it -> [it.getFileName().toString() - '.R1_val_1.fq.gz', it]  }
  .set { trimpe_r1_tuples}
trimpe_r2s
  .map { it -> [it.getFileName().toString() - '.R2_val_2.fq.gz', it]  }
  .set { trimpe_r2_tuples }

trimpe_r1_tuples
  .mix(trimpe_r2_tuples)
  .groupTuple(by: 0, sort:true )
  .map { it -> it.flatten() }
  .set { trimread }

process align {
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'

  input:
  set pair_id, file(fq1), file(fq2) from trimread
  output:
  set pair_id, file("${pair_id}.bam") into aligned
  set pair_id, file("${pair_id}.bam") into aligned2
  """
  bash $baseDir/process_scripts/alignment/dnaseqalign.sh -r $index_path -p $pair_id -x $fq1 -y $fq2 $alignopts
  """
 }

process markdups_consensus {
  publishDir "$params.output/$pair_id", mode: 'copy'

  input:
  set pair_id, file(sbam) from aligned
  output:
  set pair_id, file("${pair_id}.consensus.bam") into deduped
  script:
  """
  bash $baseDir/process_scripts/alignment/markdups.sh -a fgbio_umi -b $sbam -p $pair_id
  mv ${pair_id}.dedup.bam ${pair_id}.consensus.bam
  """
}
process markdups_picard {
  //publishDir "$params.output/$pair_id", mode: 'copy'

  input:
  set pair_id, file(sbam) from aligned2
  output:
  set pair_id, file("${pair_id}.dedup.bam") into qcbam
  script:
  """
  bash $baseDir/process_scripts/alignment/markdups.sh -a picard_umi -b $sbam -p $pair_id
  """
}

process seqqc {
  errorStrategy 'ignore'
  publishDir "$params.output/$pair_id", mode: 'copy'

  input:
  set pair_id, file(sbam) from qcbam
  output:
  file("${pair_id}.flagstat.txt") into alignstats
  file("${pair_id}.ontarget.flagstat.txt") into ontarget
  file("${pair_id}.dedupcov.txt") into dedupcov
  file("${pair_id}.meanmap.txt") into meanmap
  file("${pair_id}.libcomplex.txt") into libcomplex
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
  perl $baseDir/scripts/sequenceqc_alignment_withumi.pl -r ${index_path} *.genomecov.txt
  perl $baseDir/scripts/covstats.pl *.mapqualcov.txt
  Rscript $baseDir/scripts/plot_hist_genocov.R
  """
}

process gatkbam {
  //errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'

  input:
  set pair_id, file(sbam) from deduped

  output:
  set pair_id,file("${pair_id}.final.bam"),file("${pair_id}.final.bai") into gatkbam
  
  script:
  """
  bash $baseDir/process_scripts/variants/gatkrunner.sh -a gatkbam -b $sbam -r ${index_path} -p $pair_id
  """
}
