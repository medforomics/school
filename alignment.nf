#!/usr/bin/env nextflow

params.input = './fastq'
params.output = './analysis'

params.fastqs="$params.input/*.fastq.gz"
params.design="$params.input/design.txt"

params.genome="/project/shared/bicf_workflow_ref/GRCh38"
params.capture="$params.genome/UTSWV2.bed"
params.pairs="pe"
params.cancer="detect"

reffa=file("$params.genome/genome.fa")
dbsnp="$params.genome/dbSnp.vcf.gz"
indel="$params.genome/GoldIndels.vcf.gz"

fastqs=file(params.fastqs)
design_file = file(params.design)
dbsnp=file(dbsnp)
knownindel=file(indel)
index_path = file(params.genome)
capture_bed = file(params.capture)

snpeff_vers = 'GRCh38.82';
if (params.genome == '/project/shared/bicf_workflow_ref/GRCm38') {
   snpeff_vers = 'GRCm38.82';
}
if (params.genome == '/project/shared/bicf_workflow_ref/GRCh37') {
   snpeff_vers = 'GRCh37.75';
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
    mrgidx = header.findIndexOf{it == 'SampleMergeName'};
    oneidx = header.findIndexOf{it == 'FullPathToFqR1'};
    twoidx = header.findIndexOf{it == 'FullPathToFqR2'};
    if (twoidx == -1) {
       twoidx = oneidx
       }
    while (line = reader.readLine()) {
    	   def row = line.split("\t")
    if (fileMap.get(row[oneidx]) != null) {
	prefix << row[prefixidx]
        mername << row[mrgidx]
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
  .from(mername)
  .set { merg }

  
Channel
  .from(1..10)
  .set {counter}

Channel
  .from(1..100)
  .set {adid}


process trimpe {
  errorStrategy 'ignore'
  //publishDir "$params.output", mode: 'copy'
  input:
  val pair_ids from read_pe.buffer(size: 24, remainder: true)
  val mnames from merg.buffer(size: 24, remainder: true)
  file(read1s) from read1.buffer(size: 24, remainder: true)
  file(read2s) from read2.buffer(size: 24, remainder: true)
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
    cmd +="mv ${read1s[i]} ${mnames[i]}.batch${ct}_${i}.R1.fq.gz\n"
    cmd +="mv ${read2s[i]} ${mnames[i]}.batch${ct}_${i}.R2.fq.gz\n"
    cmd +="trim_galore --paired --stringency 3 -q 25 --illumina --gzip --length 35 ${mnames[i]}.batch${ct}_${i}.R1.fq.gz ${mnames[i]}.batch${ct}_${i}.R2.fq.gz &\n"
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
  .set { trimpe }

process alignpe {
  errorStrategy 'ignore'
  //publishDir "$params.output", mode: 'copy'

  input:
  set pair_id, file(fq1), file(fq2) from trimpe
  output:
  set val("${fq1.baseName.split("\\.batch", 2)[0]}"), file("${pair_id}.bam") into aligned
  file("${pair_id}.libcomplex.txt") into libcomplex
  when:
  params.pairs == 'pe'
  script:
  """
  source /etc/profile.d/modules.sh
  module load bwa/intel/0.7.15 samtools/intel/1.3 speedseq/20160506 picard/1.131
  bwa mem -M -t \$SLURM_CPUS_ON_NODE -R '@RG\\tID:${fq1.baseName.split("\\.batch", 2)[0]}\\tLB:tx\\tPL:illumina\\tPU:barcode\\tSM:${fq1.baseName.split("\\.batch", 2)[0]}' ${reffa} ${fq1} ${fq2} | samtools view -1 - > output.unsort.bam
  sambamba sort --tmpdir=./ -t \$SLURM_CPUS_ON_NODE -o output.sort.bam output.unsort.bam
  java -Djava.io.tmpdir=./ -Xmx20g -jar \$PICARD/picard.jar MarkDuplicatesWithMateCigar I=output.sort.bam O=${pair_id}.bam M=${pair_id}.libcomplex.txt ASSUME_SORTED=true MINIMUM_DISTANCE=300
  """
 }

aligned
   .groupTuple(by:0)
   .set {bamgrp}

process mergebam {
   publishDir "$params.output", mode: 'copy'

  input:
  set pair_id,file(bams) from bamgrp
  output:
  file("${pair_id}*hla.*") into hla
  file("${pair_id}.hist.txt") into insertsize
  set pair_id,file("${pair_id}.bam"),file("${pair_id}.bam.bai") into qcbam
  set pair_id,file("${pair_id}.bam"),file("${pair_id}.bam.bai") into targetbam
  script:
  """
  source /etc/profile.d/modules.sh
  module load samtools/intel/1.3 speedseq/20160506 picard/1.131 bwa/intel/0.7.15 bwakit/0.7.15 
  which sambamba
  count=\$(ls *.bam |wc -l)
  if [ \$count -gt 1 ]
  then
  sambamba merge -t \$SLURM_CPUS_ON_NODE merge.bam *.bam
  else
  mv *.bam merge.bam
  fi
  sambamba sort --tmpdir=./ -t \$SLURM_CPUS_ON_NODE -o ${pair_id}.bam merge.bam
  sambamba sort --tmpdir=./ -N -t \$SLURM_CPUS_ON_NODE -o output.nsort.bam merge.bam
  java -Djava.io.tmpdir=./ -Xmx4g -jar \$PICARD/picard.jar CollectInsertSizeMetrics INPUT=${pair_id}.bam HISTOGRAM_FILE=${pair_id}.hist.ps REFERENCE_SEQUENCE=${reffa} OUTPUT=${pair_id}.hist.txt
  samtools view output.nsort.bam | k8 /cm/shared/apps/bwa/intel/0.7.15/bwakit/bwa-postalt.js -p ${pair_id}.hla ${index_path}/hs38DH.fa.alt &> tmp
  run-HLA ${pair_id}.hla > ${pair_id}.hla.top 2> ${pair_id}.hla.log
  touch ${pair_id}.hla.HLA-dummy.gt
  cat ${pair_id}.hla.HLA*.gt | grep ^GT | cut -f2- > ${pair_id}.hla.all
  """
}

process seqqc {
  errorStrategy 'ignore'
  //publishDir "$params.output", mode: 'copy'

  input:
  set pair_id, file(sbam),file(idx) from qcbam
  output:
  file("${pair_id}.flagstat.txt") into alignstats
  file("${pair_id}.ontarget.flagstat.txt") into ontarget
  file("${pair_id}.mapqualcov.txt") into mapqualcov
  file("${pair_id}.dedupcov.txt") into dedupcov
  file("${pair_id}.meanmap.txt") into meanmap
  file("${pair_id}.libsizeest.txt") into libsize
  file("${pair_id}.alignmentsummarymetrics.txt") into alignmentsummarymetrics
  set file("${pair_id}_fastqc.zip"),file("${pair_id}_fastqc.html") into fastqc
  set pair_id,file("${pair_id}.ontarget.bam"),file("${pair_id}.ontarget.bam.bai") into genocovbam

  script:
  """
  source /etc/profile.d/modules.sh
  module load bedtools/2.25.0 picard/1.131 samtools/intel/1.3 fastqc/0.11.2 speedseq/20160506
  fastqc -f bam ${sbam}
  sambamba flagstat -t 30 ${sbam} > ${pair_id}.flagstat.txt
  sambamba view -t 30 -f bam -L  ${capture_bed} -o ${pair_id}.ontarget.bam ${sbam}
  sambamba flagstat -t 30 ${pair_id}.ontarget.bam > ${pair_id}.ontarget.flagstat.txt
  samtools view -b -q 1 ${pair_id}.ontarget.bam | bedtools coverage -sorted -hist -g ${index_path}/genomefile.txt -b stdin -a ${capture_bed}  >  ${pair_id}.mapqualcov.txt
  samtools view -b -F 1024 ${pair_id}.ontarget.bam | bedtools coverage -sorted -g  ${index_path}/genomefile.txt -a ${capture_bed} -b stdin -hist | grep ^all > ${pair_id}.dedupcov.txt 
  java -Djava.io.tmpdir=./ -Xmx32g -jar \$PICARD/picard.jar CollectAlignmentSummaryMetrics R=${reffa} I=${pair_id}.ontarget.bam OUTPUT=${pair_id}.alignmentsummarymetrics.txt
  java -Djava.io.tmpdir=./ -Xmx32g -jar \$PICARD/picard.jar EstimateLibraryComplexity I=${pair_id}.ontarget.bam OUTPUT=${pair_id}.libsizeest.txt
  samtools view -F 1024 ${pair_id}.ontarget.bam | awk '{sum+=\$5} END { print "Mean MAPQ =",sum/NR}' > ${pair_id}.meanmap.txt
  """
}
process genocov {
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'

  input:
  set pair_id, file(sbam),file(idx) from genocovbam
  output:
  file("${pair_id}.genomecov.txt") into genomecov
  file("${pair_id}.covhist.txt") into covhist
  file("*coverage.txt") into capcovstat
  script:
  """
  source /etc/profile.d/modules.sh
  module load bedtools/2.25.0 picard/1.131 samtools/intel/1.3 fastqc/0.11.2 speedseq/20160506
  bedtools coverage -sorted -hist -g ${index_path}/genomefile.txt -b ${sbam} -a ${capture_bed} >  ${pair_id}.covhist.txt
  perl $baseDir/scripts/calculate_depthcov.pl ${pair_id}.covhist.txt
  grep ^all ${pair_id}.covhist.txt >  ${pair_id}.genomecov.txt
  """
}

process parse_stat {
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'

  input:
  file(txt) from alignstats.toList()
  file(lc) from libcomplex.toList()
  file(ls) from libsize.toList()
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
  publishDir "$params.output", mode: 'copy'

  input:
  set pair_id, file(dbam), file(idx) from targetbam

  output:
  set pair_id,file("${pair_id}.final.bam"),file("${pair_id}.final.bai") into gatkbam
  
  script:
  """
  source /etc/profile.d/modules.sh
  module load gatk/3.5 samtools/intel/1.3
  samtools index ${dbam}
  java -Djava.io.tmpdir=./ -Xmx32g -jar \$GATK_JAR -T RealignerTargetCreator -known ${knownindel} -R ${reffa} -o ${pair_id}.bam.list -I ${dbam} -nt \$SLURM_CPUS_ON_NODE -nct 1
  java -Djava.io.tmpdir=./ -Xmx32g -jar \$GATK_JAR -I ${dbam} -R ${reffa} --filter_mismatching_base_and_quals -T IndelRealigner -targetIntervals ${pair_id}.bam.list -o ${pair_id}.realigned.bam
  java -Djava.io.tmpdir=./ -Xmx32g -jar \$GATK_JAR -l INFO -R ${reffa} --knownSites ${dbsnp} -I ${pair_id}.realigned.bam -T BaseRecalibrator -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -o ${pair_id}.recal_data.grp -nt 1 -nct \$SLURM_CPUS_ON_NODE
  java -Djava.io.tmpdir=./ -Xmx32g -jar \$GATK_JAR -T PrintReads -R ${reffa} -I ${pair_id}.realigned.bam -BQSR ${pair_id}.recal_data.grp -o ${pair_id}.final.bam -nt 1 -nct 8
  """
}
