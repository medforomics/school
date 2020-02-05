#!/usr/bin/env nextflow

params.input = './fastq'
params.output = './analysis'

params.fastqs="$params.input/*.fastq.gz"
params.design="$params.input/design.txt"

params.genome="/project/shared/bicf_workflow_ref/human/grch38_cloud/dnaref"

params.pairs="pe"
params.markdups='fgbio_umi'
params.aligner='bwa'
params.cnv = "detect"

reffa=file("$params.genome/genome.fa")
dbsnp="$params.genome/dbSnp.vcf.gz"
indel="$params.genome/GoldIndels.vcf.gz"

fastqs=file(params.fastqs)
design_file=file(params.design)
dbsnp=file(dbsnp)
knownindel=file(indel)
index_path=file(params.genome)
capture_bed = file("$params.capture")
capturedir = file("$params.capturedir")

if(capturedir.isEmpty()) {
  skipCNV = true
}
else {
  skipCNV = false
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
  set subjid,pair_id, file("${pair_id}.trim.R1.fastq.gz"),file("${pair_id}.trim.R2.fastq.gz") into trimread
  file("${pair_id}.trimreport.txt") into trimstat 
  script:
  """
  bash $baseDir/process_scripts/preproc_fastq/trimgalore.sh -p ${pair_id} -a ${read1} -b ${read2}
  perl $baseDir/process_scripts/preproc_fastq/parse_trimreport.pl ${pair_id}.trimreport.txt *trimming_report.txt
  """
}

process dalign {
  queue '32GB,128GB,256GB,256GBv1'
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid/$pair_id", mode: 'copy'

  input:
  set subjid,pair_id, file(fq1), file(fq2) from trimread
  output:
  set subjid,pair_id, file("${pair_id}.bam") into aligned
  set subjid,pair_id, file("${pair_id}.bam") into aligned2
  set subjid,pair_id, file("${pair_id}.bam"), file("${pair_id}.bam.bai") into cnvbam
  file("${pair_id}.bam.bai") into baindex
  """
  bash $baseDir/process_scripts/alignment/dnaseqalign.sh -r $index_path -a $params.aligner -p $pair_id -x $fq1 -y $fq2 $alignopts
  """
 }

process qcbam           {
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid/$pair_id", mode: 'copy'
  queue '128GB,256GB,256GBv1'

  input:
  set subjid, pair_id, file(sbam) from aligned2
  output:
  file("*fastqc*") into fastqc
  file("${pair_id}.flagstat.txt") into alignstats
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
  bash $baseDir/process_scripts/alignment/bamqc.sh -c $capture_bed -n dna -r $index_path -b $sbam -p $pair_id  
  """
}

process markdups_consensus {
  errorStrategy 'ignore'
  queue '32GB,128GB,256GB,256GBv1'
  publishDir "$params.output/$subjid/$pair_id", mode: 'copy'

  input:
  set subjid, pair_id, file(sbam) from aligned
  output:
  set subjid, pair_id, file("${pair_id}.consensus.bam"),file("${pair_id}.consensus.bam.bai") into qcbam
  set subjid, pair_id, file("${pair_id}.consensus.bam"),file("${pair_id}.consensus.bam.bai") into togatkbam
  file("*.txt") into statfile
  script:
  """
  bash $baseDir/process_scripts/alignment/markdups.sh -a $params.markdups -b $sbam -p $pair_id -r $index_path
  mv ${pair_id}.dedup.bam ${pair_id}.consensus.bam
  mv ${pair_id}.dedup.bam.bai ${pair_id}.consensus.bam.bai
  """
}
process qc_consensus {
  errorStrategy 'ignore'
  queue '128GB,256GB,256GBv1'
  publishDir "$params.output/$subjid/$pair_id", mode: 'copy'

  input:
  set subjid, pair_id, file(sbam),file(idx) from qcbam
  output:
  file("${pair_id}.uniqhist.txt") into umihist
  file("${pair_id}.dedupcov.txt") into dedupcov
  file("${pair_id}.covuniqhist.txt") into covuniqhist
  file("*coverageuniq.txt") into covuniqstat
  script:
  """
  bash $baseDir/process_scripts/alignment/bamqc.sh -c $capture_bed -n dna -r $index_path -b ${sbam} -p $pair_id -s 1
  mv ${pair_id}.genomecov.txt ${pair_id}.dedupcov.txt
  mv ${pair_id}.covhist.txt ${pair_id}.covuniqhist.txt
  mv ${pair_id}.hist.txt ${pair_id}.uniqhist.txt
  mv ${pair_id}_lowcoverage.txt ${pair_id}_lowcoverageuniq.txt
  mv ${pair_id}_exoncoverage.txt ${pair_id}_exoncoverageuniq.txt
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
  file("${pair_id}.itdseek_tandemdup.pass.vcf.gz") into itdseekvcf
  when:
  skipCNV == false
  script:
  """
  source /etc/profile.d/modules.sh
  module load htslib/gcc/1.8 
  bash $baseDir/scripts/determine_pon.sh -b $sbam -r $index_path -d $capturedir -p $pair_id
  bash $baseDir/process_scripts/variants/itdseek.sh -b $sbam -r $index_path -p $pair_id -l ${index_path}/itd_genes.bed
  perl $baseDir/process_scripts/variants/filter_itdseeker.pl -t ${pair_id} -d ${pair_id}.itdseek_tandemdup.vcf.gz
  bgzip ${pair_id}.itdseek_tandemdup.pass.vcf
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

process concatstat {
  queue '32GB,128GB,256GB,256GBv1'
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'

  input:
  file(txt) from alignstats.toList()
  file(lc) from libcomplex.toList()
  file(is) from insertsize.toList()
  file(gc) from genomecov.toList()
  //file(on) from ontarget.toList()
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
  Rscript $baseDir/scripts/plot_hist_genocov.R
  """
}
