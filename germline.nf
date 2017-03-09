#!/usr/bin/env nextflow

params.input = './fastq'
params.output = './'

params.fastqs="$params.input/*.fastq.gz"
params.design="$params.input/design.txt"

params.genome="/project/shared/bicf_workflow_ref/GRCh38"
params.capture="$params.genome/UTSWV2.bed"
params.pairs="pe"
params.cancer="detect"

gatkref=file("$params.genome/genome.fa")
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
new File(params.design).withReader { reader ->
    def hline = reader.readLine()
    def header = hline.split("\t")
    prefixidx = header.findIndexOf{it == 'SampleID'};
    oneidx = header.findIndexOf{it == 'FullPathToFqR1'};
    twoidx = header.findIndexOf{it == 'FullPathToFqR2'};
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
  publishDir "$params.output", mode: 'copy'
  input:
  val pair_ids from read_pe.buffer(size: 30, remainder: true)
  file(read1s) from read1.buffer(size: 30, remainder: true)
  file(read2s) from read2.buffer(size: 30, remainder: true)
  val ct from counter
  
  output:
  file("*_val_1.fq.gz") into trimpe_r1s mode flatten
  file("*_val_2.fq.gz") into trimpe_r2s mode flatten
  file("trimreport_${ct}.txt") into trimstat mode flatten
  script:
  assert pair_ids.size() == read1s.size()
  assert pair_ids.size() == read2s.size()
  def cmd = ''
  for( int i=0; i<pair_ids.size(); i++){
    cmd +="mv ${read1s[i]} ${pair_ids[i]}.R1.fq.gz\n"
    cmd +="mv ${read2s[i]} ${pair_ids[i]}.R2.fq.gz\n"
    cmd +="trim_galore --paired --stringency 3 -q 25 --illumina --gzip --length 35 ${pair_ids[i]}.R1.fq.gz  ${pair_ids[i]}.R2.fq.gz &\n"
  }

  """
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
  publishDir "$params.output", mode: 'copy'

  input:
  set pair_id, file(fq1), file(fq2) from trimpe
  output:
  set pair_id, file("${pair_id}.bam"), file("${pair_id}.bam.bai") into aligned
  set pair_id,file("${pair_id}.bam"),file("${pair_id}.discordants.bam"),file("${pair_id}.splitters.bam") into svbam
  file("${pair_id}.hist.txt") into insertsize
  file("${pair_id}.hla.top") into hlape
  file("${pair_id}.libcomplex.txt") into libcomplex
  when:
  params.pairs == 'pe'
  script:
  """
  module load bwakit/0.7.15 seqtk/1.2-r94 samtools/intel/1.3 speedseq/20160506 picard/1.127
  seqtk mergepe ${fq1} ${fq2} | bwa mem -M -p -t \$SLURM_CPUS_ON_NODE -R '@RG\\tID:${pair_id}\\tLB:tx\\tPL:illumina\\tPU:barcode\\tSM:${pair_id}' ${gatkref} - 2> log.bwamem | k8 /cm/shared/apps/bwakit/0.7.15/bwa-postalt.js -p ${pair_id}.hla ${gatkref}.alt | samtools view -1 - > output.unsort.bam
  run-HLA ${pair_id}.hla > ${pair_id}.hla.top 2> ${pair_id}.log.hla
  touch ${pair_id}.hla.HLA-dummy.gt
  cat ${pair_id}.hla.HLA*.gt | grep ^GT | cut -f2- > ${pair_id}.hla.all
  sambamba sort -t \$SLURM_CPUS_ON_NODE -o output.dups.bam output.unsort.bam
  java -Xmx20g -jar \$PICARD/picard.jar MarkDuplicatesWithMateCigar I=output.dups.bam O=${pair_id}.bam M=${pair_id}.libcomplex.txt ASSUME_SORTED=true MINIMUM_DISTANCE=150
  sambamba index -t 20 ${pair_id}.bam
  sambamba view -h output.unsort.bam | samblaster -M -a --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 -d discordants.sam -s splitters.sam > temp.sam
  gawk '{ if (\$0~"^@") { print; next } else { \$10="*"; \$11="*"; print } }' OFS="\\t" splitters.sam | samtools  view -S -b - | samtools sort -o ${pair_id}.splitters.bam -
  gawk '{ if (\$0~"^@") { print; next } else { \$10="*"; \$11="*"; print } }' OFS="\\t" discordants.sam | samtools  view -S  -b - | samtools sort -o ${pair_id}.discordants.bam -
  java -Xmx16g -jar \$PICARD/picard.jar CollectInsertSizeMetrics INPUT=${pair_id}.bam HISTOGRAM_FILE=${pair_id}.hist.ps REFERENCE_SEQUENCE=${gatkref} OUTPUT=${pair_id}.hist.txt
  """
}

// Calculate Metrics of Quality of Alignment
process seqqc {
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'

  input:
  set pair_id, file(sbam),file(idx) from aligned
  output:
  file("${pair_id}.flagstat.txt") into alignstats
  file("${pair_id}.subset.libcomplex.txt") into subsetlibcomplex
  file("${pair_id}.ontarget.flagstat.txt") into ontarget
  set pair_id, file("${pair_id}.ontarget.bam"),file("${pair_id}.ontarget.bam.bai") into targetbam
  set pair_id, file("${pair_id}.ontarget.bam"),file("${pair_id}.ontarget.bam.bai") into ssbam
  file("${pair_id}.genomecov.txt") into genomecov
  file("${pair_id}.mapqualcov.txt") into mapqualcov
  file("${pair_id}.dedupcov.txt") into dedupcov
  file("${pair_id}.covhist.txt") into covhist
  file("${pair_id}.meanmap.txt") into meanmap
  file("${pair_id}.alignmentsummarymetrics.txt") into alignmentsummarymetrics
  set file("${pair_id}_fastqc.zip"),file("${pair_id}_fastqc.html") into fastqc

  script:
  """
  module load bedtools/2.25.0 picard/1.127 samtools/intel/1.3 fastqc/0.11.2 speedseq/20160506
  fastqc -f bam ${sbam}
  sambamba flagstat -t 30 ${sbam} > ${pair_id}.flagstat.txt
  sambamba view -t 30 -f bam -L  ${capture_bed} -o ${pair_id}.ontarget.bam ${sbam}
  sambamba flagstat -t 30 ${pair_id}.ontarget.bam > ${pair_id}.ontarget.flagstat.txt
  bedtools coverage -sorted -hist -g ${index_path}/genomefile.txt -b ${pair_id}.ontarget.bam -a ${capture_bed} | grep ^all >  ${pair_id}.covhist.txt
  grep ^all ${pair_id}.covhist.txt >  ${pair_id}.genomecov.txt
  samtools view -b -q 1 ${pair_id}.ontarget.bam | bedtools coverage -sorted -hist -g ${index_path}/genomefile.txt -b stdin -a ${capture_bed}  >  ${pair_id}.mapqualcov.txt
  samtools view -b -F 1024 ${pair_id}.ontarget.bam | bedtools coverage -sorted -g  ${index_path}/genomefile.txt -a ${capture_bed} -b stdin -hist | grep ^all > ${pair_id}.dedupcov.txt 
  perl $baseDir/scripts/subset_bam.pl ${pair_id}.ontarget.bam ${pair_id}.ontarget.flagstat.txt ${capture_bed}
  java -Xmx8g -jar \$PICARD/picard.jar EstimateLibraryComplexity INPUT=${pair_id}.subset.bam OUTPUT=${pair_id}.subset.libcomplex.txt
  java -Xmx20g -jar \$PICARD/picard.jar CollectAlignmentSummaryMetrics R=${gatkref} I=${pair_id}.ontarget.bam OUTPUT=${pair_id}.alignmentsummarymetrics.txt
  samtools view -F 1024 ${pair_id}.ontarget.bam | awk '{sum+=\$5} END { print "Mean MAPQ =",sum/NR}' > ${pair_id}.meanmap.txt
 """
}

process parse_stat {
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'

  input:
  file(txt) from alignstats.toList()
  file(lc) from libcomplex.toList()
  file(sc) from subsetlibcomplex.toList()
  file(is) from insertsize.toList()
  file(gc) from genomecov.toList()
  file(on) from ontarget.toList()
  file(tr) from trimstat.toList()
  file(mq) from mapqualcov.toList()
  file(de) from dedupcov.toList()
  file(mm) from meanmap.toList()
  file(asmet) from alignmentsummarymetrics.toList()
  
  output:
  file('sequence.stats.txt')
  file('*.png')
  script:
  """
  perl $baseDir/scripts/parse_seqqc_bulktrim.pl *.genomecov.txt
  perl $baseDir/scripts/covstats.pl *.mapqualcov.txt
  Rscript $baseDir/scripts/plot_hist_genocov.R
  """
}

process svcall {
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'
  input:
  set pair_id,file(ssbam),file(discordbam),file(splitters) from svbam
  output:
  set pair_id,file("${pair_id}.sssv.sv.vcf.gz") into svvcf
  set pair_id,file("${pair_id}.sssv.sv.vcf.gz.annot.txt") into svannot
  script:
  """
  module load samtools/intel/1.3 bedtools/2.25.0 bcftools/intel/1.3 snpeff/4.2 speedseq/20160506 vcftools/0.1.14
  speedseq sv -t 30 -o ${pair_id}.sssv -R ${gatkref} -B ${ssbam} -D ${discordbam} -S ${splitters} -x ${index_path}/exclude_alt.bed
  perl $baseDir/scripts/parse_svresults.pl -i ${pair_id}.sssv.sv.vcf.gz -r ${index_path}
  """
}

process gatkbam {
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'

  input:
  set pair_id, file(dbam), file(idx) from targetbam

  output:
  set pair_id,file("${pair_id}.final.bam"),file("${pair_id}.final.bai") into gatkbam
  set pair_id,file("${pair_id}.final.bam"),file("${pair_id}.final.bai") into sambam
  set pair_id,file("${pair_id}.final.bam"),file("${pair_id}.final.bai") into hsbam
  set pair_id,file("${pair_id}.final.bam"),file("${pair_id}.final.bai") into platbam
  
  script:
  """
  module load gatk/3.5 samtools/intel/1.3
  samtools index ${dbam}
  java -Xmx32g -jar \$GATK_JAR -T RealignerTargetCreator -known ${knownindel} -R ${gatkref} -o ${pair_id}.bam.list -I ${dbam} -nt 30 -nct 1
  java -Xmx32g -jar \$GATK_JAR -I ${dbam} -R ${gatkref} --filter_mismatching_base_and_quals -T IndelRealigner -targetIntervals ${pair_id}.bam.list -o ${pair_id}.realigned.bam
  java -Xmx32g -jar \$GATK_JAR -l INFO -R ${gatkref} --knownSites ${dbsnp} -I ${pair_id}.realigned.bam -T BaseRecalibrator -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -o ${pair_id}.recal_data.grp -nt 1 -nct 30
  java -Xmx32g -jar \$GATK_JAR -T PrintReads -R ${gatkref} -I ${pair_id}.realigned.bam -BQSR ${pair_id}.recal_data.grp -o ${pair_id}.final.bam -nt 1 -nct 8
  """
}

process gatkgvcf {
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'
  cpus 30

  input:
  set pair_id,file(gbam),file(gidx) from gatkbam
  output:
  file("${pair_id}.gatkpanel.vcf.gz") into gatkfilt
  set pair_id,file("${pair_id}.gatk.vcf.gz") into gatkvcf

  script:
  """
  module load gatk/3.5 python/2.7.x-anaconda bedtools/2.25.0 snpeff/4.2 vcftools/0.1.14 
  java -Xmx32g -jar \$GATK_JAR -R ${gatkref} -D ${dbsnp} -T HaplotypeCaller -stand_call_conf 30 -stand_emit_conf 10.0 -A FisherStrand -A QualByDepth -A VariantType -A DepthPerAlleleBySample -A HaplotypeScore -A AlleleBalance -variant_index_type LINEAR -variant_index_parameter 128000 --emitRefConfidence GVCF -I ${gbam} -o ${pair_id}.gatk.g.vcf -nct 2
  java -Xmx32g -jar \$GATK_JAR -R ${gatkref} -D ${dbsnp} -T GenotypeGVCFs -o gatk.vcf -nt 4 --variant ${pair_id}.gatk.g.vcf
  vcf-annotate -n --fill-type gatk.vcf | bcftools norm -c s -f ${gatkref} -w 10 -O z -o ${pair_id}.gatk.vcf.gz -
  tabix ${pair_id}.gatk.vcf.gz
  java -Xmx32g -jar \$GATK_JAR -R ${gatkref} -T VariantFiltration -V ${pair_id}.gatk.vcf.gz -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o ${pair_id}.filtgatk.vcf 
  java -jar \$SNPEFF_HOME/SnpSift.jar filter '((QUAL >= 10) & (MQ > 40) & (DP >= 10))' ${pair_id}.filtgatk.vcf | bedtools intersect -header -a stdin -b ${capture_bed} |bgzip > ${pair_id}.gatkpanel.vcf.gz
  """
}

process mpileup {
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'
  cpus 32

  input:
  set pair_id,file(gbam),file(gidx) from sambam
  output:
  file("${pair_id}.sampanel.vcf.gz") into samfilt
  set pair_id,file("${pair_id}.sam.vcf.gz") into samvcf
  script:
  """
  module load python/2.7.x-anaconda samtools/intel/1.3 bedtools/2.25.0 bcftools/intel/1.3 snpeff/4.2 vcftools/0.1.14
  samtools mpileup -t 'AD,DP,INFO/AD' -ug -Q20 -C50 -f ${gatkref} ${gbam} | bcftools call -vmO z -o ${pair_id}.sam.ori.vcf.gz
  vcf-concat ${pair_id}.sam.ori.vcf.gz | vcf-sort |vcf-annotate -n --fill-type | bcftools norm -c s -f ${gatkref} -w 10 -O z -o ${pair_id}.sam.vcf.gz -
  java -jar \$SNPEFF_HOME/SnpSift.jar filter '((QUAL >= 10) & (MQ >= 20) & (DP >= 10))' ${pair_id}.sam.vcf.gz |bedtools intersect -header -a stdin -b ${capture_bed} |bgzip > ${pair_id}.sampanel.vcf.gz
  """
}
process hotspot {
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'

  input:
  set pair_id,file(gbam),file(gidx) from hsbam
  output:
  set pair_id,file("${pair_id}.hotspot.vcf.gz") into hsvcf
  when:
  params.cancer == "detect"
  script:
  """
  module load python/2.7.x-anaconda samtools/intel/1.3 bedtools/2.25.0 bcftools/intel/1.3 snpeff/4.2 vcftools/0.1.14
  samtools mpileup -d 99999 -t 'AD,DP,INFO/AD' -uf ${gatkref} ${gbam} > ${pair_id}.mpi
  bcftools filter -i "AD[1]/DP > 0.01" ${pair_id}.mpi | bcftools filter -i "DP > 50" | bcftools call -m -A |vcf-annotate -n --fill-type |  bcftools norm -c s -f /project/shared/bicf_workflow_ref/GRCh38/genome.fa -w 10 -O z -o ${pair_id}.lowfreq.vcf.gz -
  java -jar \$SNPEFF_HOME/SnpSift.jar annotate ${index_path}/cosmic.vcf.gz ${pair_id}.lowfreq.vcf.gz | java -jar \$SNPEFF_HOME/SnpSift.jar filter "(CNT[*] >0)" - |bgzip > ${pair_id}.hotspot.vcf.gz
  """
}
process speedseq {
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'

  input:
  set pair_id,file(gbam),file(gidx) from ssbam
  output:
  file("${pair_id}.sspanel.vcf.gz") into ssfilt
  set pair_id,file("${pair_id}.ssvar.vcf.gz") into ssvcf

  script:
  """
  module load python/2.7.x-anaconda samtools/intel/1.3 bedtools/2.25.0 bcftools/intel/1.3 snpeff/4.2 speedseq/20160506 vcftools/0.1.14
  speedseq var -t \$SLURM_CPUS_ON_NODE -o ssvar ${gatkref} ${gbam}
  vcf-annotate -n --fill-type ssvar.vcf.gz| bcftools norm -c s -f ${gatkref} -w 10 -O z -o ${pair_id}.ssvar.vcf.gz -
  java -Xmx10g -jar \$SNPEFF_HOME/SnpSift.jar filter '((QUAL >= 10) & (DP >= 10))' ${pair_id}.ssvar.vcf.gz |bedtools intersect -header -a stdin -b ${capture_bed} | bgzip > ${pair_id}.sspanel.vcf.gz
  """
}
process platypus {
  errorStrategy 'ignore'

  publishDir "$params.output", mode: 'copy'
  cpus 32

  input:
  set pair_id,file(gbam),file(gidx) from platbam

  output:
  file("${pair_id}.platpanel.vcf.gz") into platfilt
  set pair_id,file("${pair_id}.platypus.vcf.gz") into platvcf

  script:
  """
  module load python/2.7.x-anaconda bedtools/2.25.0 snpeff/4.2 platypus/gcc/0.8.1 bcftools/intel/1.3 samtools/intel/1.3 vcftools/0.1.14
  Platypus.py callVariants --minMapQual=10 --mergeClusteredVariants=1 --nCPU=30 --bamFiles=${gbam} --refFile=${gatkref} --output=platypus.vcf
  bgzip platypus.vcf
  tabix platypus.vcf.gz
  vcf-sort platypus.vcf.gz| vcf-annotate -n --fill-type -n | bcftools norm -c s -f ${gatkref} -w 10 -O z -o ${pair_id}.platypus.vcf.gz -
  java -jar \$SNPEFF_HOME/SnpSift.jar filter "((QUAL >= 10) & (QD > 2) & (FILTER = 'PASS'))" ${pair_id}.platypus.vcf.gz | bedtools intersect -header -a stdin -b ${capture_bed} |bgzip > ${pair_id}.platpanel.vcf.gz
  """
}

ssvcf .phase(gatkvcf)
      .map {p,q -> [p[0],p[1],q[1]]}
      .set { twovcf }
twovcf .phase(samvcf)
      .map {p,q -> [p[0],p[1],p[2],q[1]]}
      .set { threevcf }
threevcf .phase(platvcf)
      .map {p,q -> [p[0],p[1],p[2],p[3],q[1]]}
      .set { fourvcf }
if (params.cancer == "detect") {
  fourvcf .phase(hsvcf)
  	.map {p,q -> [p[0],p[1],p[2],p[3],p[4],q[1]]}
      	.set { vcflist }
}
else {
  Channel
	.from(fourvcf)
  	.into {vcflist}
}

process integrate {
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'

  input:
  set fname,file(ss),file(gatk),file(sam),file(plat),file(hs) from vcflist
  
  output:
  set fname,file("${fname}.union.vcf.gz") into union
  script:
  if (params.cancer == "detect")
  """
  module load gatk/3.5 python/2.7.x-anaconda bedtools/2.25.0 snpeff/4.2 bcftools/intel/1.3 samtools/intel/1.3 vcftools/0.1.14
  bedtools multiinter -i ${gatk} ${sam} ${ss} ${plat} ${hs} -names gatk sam ssvar platypus hotspot |cut -f 1,2,3,5 | bedtools sort -i stdin | bedtools merge -c 4 -o distinct >  ${fname}_integrate.bed
  bedtools intersect -header -v -a ${hs} -b ${sam} |bedtools intersect -header -v -a stdin -b ${gatk} | bedtools intersect -header -v -a stdin -b ${ss} |  bedtools intersect -header -v -a stdin -b ${plat} | bgzip > ${fname}.hotspot.nooverlap.vcf.gz
  tabix ${fname}.hotspot.nooverlap.vcf.gz
  tabix ${gatk}
  tabix ${sam}
  tabix ${ss}
  tabix ${plat}
  java -Xmx32g -jar \$GATK_JAR -R ${gatkref} -T CombineVariants --filteredrecordsmergetype KEEP_UNCONDITIONAL --variant:gatk ${gatk} --variant:sam ${sam} --variant:ssvar ${ss} --variant:platypus ${plat} --variant:hotspot ${fname}.hotspot.nooverlap.vcf.gz -genotypeMergeOptions PRIORITIZE -priority sam,ssvar,gatk,platypus,hotspot -o ${fname}.int.vcf
  perl $baseDir/scripts/uniform_integrated_vcf.pl ${fname}.int.vcf
  bgzip ${fname}_integrate.bed
  tabix ${fname}_integrate.bed.gz
  bgzip ${fname}.uniform.vcf
  tabix ${fname}.uniform.vcf.gz
  bcftools annotate -a ${fname}_integrate.bed.gz --columns CHROM,FROM,TO,CallSet -h ${index_path}/CallSet.header ${fname}.uniform.vcf.gz | bgzip > ${fname}.union.vcf.gz
  """
  else
  """
  module load gatk/3.5 python/2.7.x-anaconda bedtools/2.25.0 snpeff/4.2 bcftools/intel/1.3 samtools/intel/1.3
  module load vcftools/0.1.14
  bedtools multiinter -i ${gatk} ${sam} ${ss} ${plat} -names gatk sam ssvar platypus |cut -f 1,2,3,5 | bedtools sort -i stdin | bedtools merge -c 4 -o distinct >  ${fname}_integrate.bed
  java -Xmx32g -jar \$GATK_JAR -R ${gatkref} -T CombineVariants --filteredrecordsmergetype KEEP_UNCONDITIONAL -genotypeMergeOptions PRIORITIZE --variant:gatk ${gatk} --variant:sam ${sam} --variant:ssvar ${ss} --variant:platypus ${plat} -priority sam,ssvar,gatk,platypus -o ${fname}.int.vcf
  perl $baseDir/scripts/uniform_integrated_vcf.pl ${fname}.int.vcf
  bgzip ${fname}_integrate.bed
  tabix ${fname}_integrate.bed.gz
  bgzip ${fname}.uniform.vcf
  tabix ${fname}.uniform.vcf.gz
  bcftools annotate -a ${fname}_integrate.bed.gz --columns CHROM,FROM,TO,CallSet -h ${index_path}/CallSet.header ${fname}.uniform.vcf.gz | bgzip > ${fname}.union.vcf.gz
  """
}

process annot {
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'

  input:
  set fname,unionvcf from union
  
  output:
  file("${fname}.annot.vcf.gz") into annotvcf
  file("${fname}.annot.vcf.gz") into coding
  file("${fname}.stats.txt") into stats
  file("${fname}.statplot*") into plotstats

  script:
  """
  module load python/2.7.x-anaconda bedtools/2.25.0 snpeff/4.2 bcftools/intel/1.3 samtools/intel/1.3
  tabix ${unionvcf}
  bcftools annotate -Oz -a ${index_path}/ExAC.vcf.gz -o ${fname}.exac.vcf.gz --columns CHROM,POS,AC_Het,AC_Hom,AC_Hemi,AC_Adj,AN_Adj,AC_POPMAX,AN_POPMAX,POPMAX ${unionvcf}
  tabix ${fname}.exac.vcf.gz 
  bcftools annotate -Oz -a ${index_path}/dbSnp.vcf.gz -o ${fname}.dbsnp.vcf.gz --columns CHROM,POS,ID,RS ${fname}.exac.vcf.gz
  tabix ${fname}.dbsnp.vcf.gz
  bcftools annotate -Oz -a ${index_path}/clinvar.vcf.gz -o ${fname}.clinvar.vcf.gz --columns CHROM,POS,CLNSIG,CLNDSDB,CLNDSDBID,CLNDBN,CLNREVSTAT,CLNACC ${fname}.dbsnp.vcf.gz
  tabix ${fname}.clinvar.vcf.gz
  java -Xmx10g -jar \$SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c \$SNPEFF_HOME/snpEff.config ${snpeff_vers} ${fname}.clinvar.vcf.gz | java -jar \$SNPEFF_HOME/SnpSift.jar annotate ${index_path}/cosmic.vcf.gz - | java -Xmx10g -jar \$SNPEFF_HOME/SnpSift.jar dbnsfp -v -db ${index_path}/dbNSFP.txt.gz - | java -Xmx10g -jar \$SNPEFF_HOME/SnpSift.jar gwasCat -db ${index_path}/gwas_catalog.tsv - |bgzip > ${fname}.annot.vcf.gz
  tabix ${fname}.annot.vcf.gz
  bcftools stats ${fname}.annot.vcf.gz > ${fname}.stats.txt
  plot-vcfstats -s -p ${fname}.statplot ${fname}.stats.txt
  zcat ${fname}.annot.vcf.gz | \$SNPEFF_HOME/scripts/vcfEffOnePerLine.pl | java -jar \$SNPEFF_HOME/SnpSift.jar extractFields - CHROM POS ID ANN[*].GENE ANN[*].HGVS_P GEN[0].AD GEN[0].DP MQ CallSet AC_POPMAX AN_POPMAX ANN[*].EFFECT ANN[*].IMPACT CNT[*]| uniq |grep -v "MODIFIER" |uniq > ${fname}.coding.txt
  java -jar \$SNPEFF_HOME/SnpSift.jar filter "(GEN[0].DP > 49 & (ANN[*].IMPACT = 'HIGH' | ANN[*].IMPACT = 'MODERATE'))" ${fname}.annot.vcf.gz | bgzip > ${fname}.coding.vcf.gz
  """
}
