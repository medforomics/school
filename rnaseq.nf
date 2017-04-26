#!/usr/bin/env nextflow
	
// Default parameter values to run tests
params.input = './fastq'
params.output = './'

params.fastqs="$params.input/*.fastq.gz"
params.design="$params.input/design.txt"

params.genome="/project/shared/bicf_workflow_ref/GRCh38/"
params.markdups="skip"
params.stranded="0"
params.pairs="pe"
params.geneset = 'h.all.v5.1.symbols.gmt'
params.align = 'hisat'
params.fusion = 'detect'
params.dea = 'skip'
params.cancer="detect"
params.variant = "detect"

design_file = file(params.design)
fastqs=file(params.fastqs)
design_file = file(params.design)
gtf_file = file("$params.genome/gencode.gtf")
stringtie_gtf=file("$params.genome/gencode.chrrm.gtf")
genenames = file("$params.genome/genenames.txt")
geneset = file("$params.genome/gsea_gmt/$params.geneset")
dbsnp="$params.genome/dbSnp.vcf.gz"
indel="$params.genome/GoldIndels.vcf.gz"
knownindel=file(indel)
dbsnp=file(dbsnp)

snpeff_vers = 'GRCh38.82';
if (params.genome == '/project/shared/bicf_workflow_ref/GRCm38') {
   snpeff_vers = 'GRCm38.82';
}
if (params.genome == '/project/shared/bicf_workflow_ref/GRCh37') {
   snpeff_vers = 'GRCh37.75';
}

// params genome is the directory
// base name for the index is always genome
index_path = file(params.genome)
index_name = "genome"
star_index = 'star_index/'

// Pair handling, helper function taken from rnatoy
// which is covered by the GNU General Public License v3
// https://github.com/nextflow-io/rnatoy/blob/master/main.nf

def fileMap = [:]

fastqs.each {
    final fileName = it.getFileName().toString()
    prefix = fileName.lastIndexOf('/')
    fileMap[fileName] = it
}
def prefix = []
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
	      prefix << tuple(row[prefixidx],fileMap.get(row[oneidx]),fileMap.get(row[twoidx]))
	   }
	  
} 
}
if( ! prefix) { error "Didn't match any input files with entries in the design file" }

if (params.pairs == 'pe') {
Channel
  .from(prefix)
  .set { read_pe }
Channel
  .empty()
  .set { read_se } 
}
if (params.pairs == 'se') {
Channel
  .from(prefix)
  .into { read_se }
Channel
  .empty()
  .set { read_pe }
}

//
// Trim raw reads using trimgalore
//
process trimpe {
  errorStrategy 'ignore'
  input:
  set pair_id, file(read1), file(read2) from read_pe
  output:
  set pair_id, file("${read1.baseName.split("\\.fastq", 2)[0]}_val_1.fq.gz"), file("${read2.baseName.split("\\.fastq", 2)[0]}_val_2.fq.gz") into trimpe
  set pair_id, file("${read1.baseName.split("\\.fastq", 2)[0]}_val_1.fq.gz"), file("${read2.baseName.split("\\.fastq", 2)[0]}_val_2.fq.gz") into fusionfq
  script:
  """
  module load trimgalore/0.4.1 cutadapt/1.9.1
  trim_galore --paired -q 25 --illumina --gzip --length 35 ${read1} ${read2}
  """
}
process trimse {
  errorStrategy 'ignore'
  input:
  set pair_id, file(read1) from read_se
  output:
  set pair_id, file("${read1.baseName.split("\\.fastq", 2)[0]}_trimmed.fq.gz") into trimse
  script:
  """
  module load trimgalore/0.4.1 cutadapt/1.9.1
  trim_galore -q 25 --illumina --gzip --length 35 ${read1}
  """
}

//
// Align trimmed reads to genome indes with hisat2
// Sort and index with samtools
// QC aligned reads with fastqc
// Alignment stats with samtools
//

process starfusion {
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'
  input:
  set pair_id, file(fq1), file(fq2) from fusionfq
  output:
  file("${pair_id}.starfusion.txt") into fusionout
  when:
  params.fusion == 'detect'
  script:
  """
  module add python/2.7.x-anaconda star/2.5.2b
  STAR-Fusion --genome_lib_dir ${index_path}/CTAT_lib/ --left_fq ${fq1} --right_fq ${fq2} --output_dir star_fusion &> star_fusion.err
  mv star_fusion/star-fusion.fusion_candidates.final.abridged ${pair_id}.starfusion.txt
  """
}

process alignpe {
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'
  cpus 32

  input:
  set pair_id, file(fq1), file(fq2) from trimpe
  output:
  set pair_id, file("${pair_id}.bam") into alignpe
  file("${pair_id}.flagstat.txt") into alignstats_pe
  file("${pair_id}.hisatout.txt") into hsatoutpe
  set file("${pair_id}_fastqc.zip"),file("${pair_id}_fastqc.html") into fastqcpe
  script:
  if (params.align == 'hisat')
  """
  module load hisat2/2.0.1-beta-intel samtools/intel/1.3 fastqc/0.11.2 picard/1.127 speedseq/20160506
  hisat2 -p 30 --no-unal --dta -x ${index_path}/${index_name} -1 ${fq1} -2 ${fq2} -S out.sam 2> ${pair_id}.hisatout.txt
  sambamba view -t 30 -f bam -S -o output.bam out.sam
  sambamba sort -t 30 -o ${pair_id}.bam output.bam
  sambamba flagstat -t 30 ${pair_id}.bam > ${pair_id}.flagstat.txt
  fastqc -f bam ${pair_id}.bam
  """
  else
  """
  module load star/2.4.2a samtools/intel/1.3 fastqc/0.11.2 picard/1.127 speedseq/20160506
  STAR --genomeDir ${index_path}/${star_index} --readFilesIn ${fq1} ${fq2} --readFilesCommand zcat --genomeLoad NoSharedMemory --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMheaderCommentFile COfile.txt --outSAMheaderHD @HD VN:1.4 SO:coordinate --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --sjdbScore 1 --limitBAMsortRAM 60000000000 --outFileNamePrefix out
  mv outLog.final.out ${pair_id}.hisatout.txt
  sambamba sort -t 30 -o ${pair_id}.bam outAligned.sortedByCoord.out.bam
  sambamba flagstat -t 30 ${pair_id}.bam > ${pair_id}.flagstat.txt
  fastqc -f bam ${pair_id}.bam
  """
}
process alignse {
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'
  cpus 32

  input:
  set pair_id, file(fq1) from trimse
  output:
  set pair_id, file("${pair_id}.bam") into alignse
  file("${pair_id}.flagstat.txt") into alignstats_se
  file("${pair_id}.hisatout.txt") into hsatoutse
  set file("${pair_id}_fastqc.zip"),file("${pair_id}_fastqc.html") into fastqcse
  script:
  if (params.align == 'hisat')
  """
  module load hisat2/2.0.1-beta-intel samtools/intel/1.3 fastqc/0.11.2 speedseq/20160506 picard/1.127
  hisat2 -p 30 --no-unal --dta -x ${index_path}/${index_name} -U ${fq1} -S out.sam 2> ${pair_id}.hisatout.txt
  sambamba view -t 30 -f bam -S -o output.bam out.sam
  sambamba sort -t 30 -o ${pair_id}.bam output.bam
  sambamba flagstat -t 30 ${pair_id}.bam > ${pair_id}.flagstat.txt
  fastqc -f bam ${pair_id}.bam
  """
  else
  """
  module load star/2.4.2a samtools/intel/1.3 fastqc/0.11.2 picard/1.127 speedseq/20160506
  STAR --genomeDir ${index_path}/${star_index} --readFilesIn ${fq1} --readFilesCommand zcat --genomeLoad NoSharedMemory --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMheaderCommentFile COfile.txt --outSAMheaderHD @HD VN:1.4 SO:coordinate --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD --outSAMstrandField intronMotif --outSAMtype SAM --quantMode TranscriptomeSAM --sjdbScore 1 --limitBAMsortRAM 60000000000 --outFileNamePrefix out
  mv outLog.final.out ${pair_id}.hisatout.txt
  sambamba view -t 30 -f bam -S -o output.bam outAligned.out.sam
  sambamba sort -t 30 -o ${pair_id}.bam output.bam
  sambamba flagstat -t 30 ${pair_id}.bam > ${pair_id}.flagstat.txt
  fastqc -f bam ${pair_id}.bam
  """
}

// From here on we are the same for PE and SE, so merge channels and carry on

Channel
  .empty()
  .mix(alignse, alignpe)
  .tap { aligned2 }
  .set { aligned }

Channel
  .empty()
  .mix(alignstats_se, alignstats_pe)
  .set { alignstats }

Channel
  .empty()
  .mix(hsatoutse, hsatoutpe)
  .set { hsatout }

//
// Summarize all flagstat output
//
process parse_alignstat {

  publishDir "$params.output", mode: 'copy'

  input:
  file(txt) from alignstats.toList()
  file(txt) from  hsatout.toList()

  output:
  file('alignment.summary.txt')

  """
  perl $baseDir/scripts/parse_flagstat.pl *.flagstat.txt
  """
}

//
// Identify duplicate reads with Picard
//
process markdups {
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'

  input:
  set pair_id, file(sbam) from aligned
  output:
  set pair_id, file("${pair_id}.dedup.bam") into deduped1
  set pair_id, file("${pair_id}.dedup.bam") into deduped2
  set pair_id, file("${pair_id}.dedup.bam"), file("${pair_id}.dedup.bam.bai") into deduped3
  script:
  if (params.markdups == 'mark')
  """
  module load picard/1.127 speedseq/20160506
  sambamba markdup -t 20 -r ${sbam} ${pair_id}.dedup.bam
  """
  else
  """
  module load samtools/intel/1.3
  cp ${sbam} ${pair_id}.dedup.bam
  samtools index ${pair_id}.dedup.bam
  """
}

process gatkbam {
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'

  input:
  set pair_id, file(rbam) from deduped3
  output:
  set pair_id,file("${pair_id}.final.bam"),file("${pair_id}.final.bai") into gatkbam
  set pair_id,file("${pair_id}.final.bam"),file("${pair_id}.final.bai") into sambam
  set pair_id,file("${pair_id}.final.bam"),file("${pair_id}.final.bai") into hsbam
  set pair_id,file("${pair_id}.final.bam"),file("${pair_id}.final.bai") into ssbam
  when:
  params.variant == "detect"
  script:
  """
  module load gatk/3.5 samtools/intel/1.3 speedseq/20160506 picard/1.127
  java -Xmx4g -jar \$PICARD/picard.jar CleanSam INPUT=${rbam} O=${pair_id}.clean.bam
  java -Xmx4g -jar \$PICARD/picard.jar AddOrReplaceReadGroups INPUT=${pair_id}.clean.bam O=${pair_id}.rg_added_sorted.bam SO=coordinate RGID=${pair_id} RGLB=tx RGPL=illumina RGPU=barcode RGSM=${pair_id}
  sambamba index ${pair_id}.rg_added_sorted.bam
  java -Xmx4g -jar \$GATK_JAR -T SplitNCigarReads -R ${index_path}/hisat_genome.fa -I ${pair_id}.rg_added_sorted.bam -o ${pair_id}.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
  java -Xmx4g -jar \$GATK_JAR -l INFO -R ${index_path}/hisat_genome.fa --knownSites ${dbsnp} -I ${pair_id}.split.bam -T BaseRecalibrator -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -o ${pair_id}.recal_data.grp -nt 1 -nct 30
  java -Xmx4g -jar \$GATK_JAR -T PrintReads -R ${index_path}/hisat_genome.fa -I ${pair_id}.split.bam -BQSR ${pair_id}.recal_data.grp -o ${pair_id}.final.bam -nt 1 -nct 8
  """
}
process mpileup {
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'
  cpus 32

  input:
  set pair_id,file(gbam),file(gidx) from sambam
  output:
  set pair_id,file("${pair_id}.sam.vcf.gz") into samvcf
  when:
  params.variant == "detect"
  script:
  """
  module load python/2.7.x-anaconda samtools/intel/1.3 bedtools/2.25.0 bcftools/intel/1.3 snpeff/4.2 vcftools/0.1.14
  which samtools
  samtools mpileup -t 'AD,DP,INFO/AD' -ug -Q20 -C50 -f ${index_path}/hisat_genome.fa ${gbam} | bcftools call -vmO z -o ${pair_id}.sam.ori.vcf.gz
  vcf-concat ${pair_id}.sam.ori.vcf.gz | vcf-sort |vcf-annotate -n --fill-type | bcftools norm -c s -f ${index_path}/hisat_genome.fa -w 10 -O z -o ${pair_id}.sam.vcf.gz -
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
  params.cancer == "detect" && params.variant == "detect"
  script:
  """
  module load python/2.7.x-anaconda samtools/intel/1.3 bedtools/2.25.0 bcftools/intel/1.3 snpeff/4.2 vcftools/0.1.14
  which samtools
  samtools mpileup -d 99999 -t 'AD,DP,INFO/AD' -uf ${index_path}/hisat_genome.fa ${gbam} > ${pair_id}.mpi
  bcftools filter -i "AD[1]/DP > 0.01" ${pair_id}.mpi | bcftools filter -i "DP > 50" | bcftools call -m -A |vcf-annotate -n --fill-type |  bcftools norm -c s -f /project/shared/bicf_workflow_ref/GRCh38/hisat_genome.fa -w 10 -O z -o ${pair_id}.lowfreq.vcf.gz -
  java -jar \$SNPEFF_HOME/SnpSift.jar annotate ${index_path}/cosmic.vcf.gz ${pair_id}.lowfreq.vcf.gz | java -jar \$SNPEFF_HOME/SnpSift.jar filter "(CNT[*] >0)" - |bgzip > ${pair_id}.hotspot.vcf.gz
  """
}
process speedseq {
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'

  input:
  set pair_id,file(gbam),file(gidx) from ssbam
  output:
  set pair_id,file("${pair_id}.ssvar.vcf.gz") into ssvcf
  when:
  params.variant == "detect"
  script:
  """
  module load python/2.7.x-anaconda samtools/intel/1.3 bedtools/2.25.0 bcftools/intel/1.3 snpeff/4.2 speedseq/20160506 vcftools/0.1.14
  speedseq var -t \$SLURM_CPUS_ON_NODE -o ssvar ${index_path}/hisat_genome.fa ${gbam}
  tabix -f ssvar.vcf.gz
  bcftools norm -c s -f /project/shared/bicf_workflow_ref/GRCh38/hisat_genome.fa -w 10 -O v ssvar.vcf.gz 2> bcftools.err | vcf-annotate -n --fill-type |bgzip > ${pair_id}.ssvar.vcf.gz
  """
}

process gatkgvcf {
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'
  cpus 30

  input:
  set pair_id,file(gbam),file(gidx) from gatkbam
  output:
  set pair_id, file("${pair_id}.gatk.vcf.gz") into gatkvcf
  when:
  params.variant == "detect"
  script:
  """
  module load python/2.7.x-anaconda gatk/3.5 bedtools/2.25.0 snpeff/4.2 vcftools/0.1.14 
  java -Xmx64g -jar \$GATK_JAR -R ${index_path}/hisat_genome.fa -D ${dbsnp} -T HaplotypeCaller -stand_call_conf 30 -stand_emit_conf 10.0 -A FisherStrand -A QualByDepth -A VariantType -A DepthPerAlleleBySample -A HaplotypeScore -A AlleleBalance -variant_index_type LINEAR -variant_index_parameter 128000 --emitRefConfidence GVCF -I ${gbam} -o ${pair_id}.gatk.g.vcf -nct 2
  java -Xmx64g -jar \$GATK_JAR -R ${index_path}/hisat_genome.fa -D ${dbsnp} -T GenotypeGVCFs -o gatk.vcf -nt 4 --variant ${pair_id}.gatk.g.vcf
  vcf-annotate -n --fill-type gatk.vcf | bcftools norm -c s -f ${index_path}/hisat_genome.fa -w 10 -O z -o ${pair_id}.gatk.vcf.gz - &> bcftools.out
  tabix ${pair_id}.gatk.vcf.gz
   """
}

if (params.variant == "detect") {
   ssvcf .phase(gatkvcf)
      .map {p,q -> [p[0],p[1],q[1]]}
      .set { twovcf }
   twovcf .phase(samvcf)
      .map {p,q -> [p[0],p[1],p[2],q[1]]}
      .set { threevcf }
   if (params.cancer == "detect") {
      threevcf .phase(hsvcf)
  	.map {p,q -> [p[0],p[1],p[2],p[3],q[1]]}
      	.set { vcflist }
   }
   else {
   	Channel
		.from(threevcf)
  		.into {vcflist}
	}
}

process integrate {
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'

  input:
  set fname,file(ss),file(gatk),file(sam),file(hs) from vcflist
  
  output:
  set fname,file("${fname}.union.vcf.gz") into union
  when:
  params.variant == "detect"
  script:
  if (params.cancer == "detect")
  """
  module load gatk/3.5 python/2.7.x-anaconda bedtools/2.25.0 snpeff/4.2 bcftools/intel/1.3 samtools/intel/1.3 vcftools/0.1.14
  bedtools multiinter -i ${gatk} ${sam} ${ss} ${hs} -names gatk sam ssvar hotspot |cut -f 1,2,3,5 | bedtools sort -i stdin | bedtools merge -c 4 -o distinct >  ${fname}_integrate.bed
  bedtools intersect -header -v -a ${hs} -b ${sam} |bedtools intersect -header -v -a stdin -b ${gatk} | bedtools intersect -header -v -a stdin -b ${ss} | bgzip > ${fname}.hotspot.nooverlap.vcf.gz
  tabix ${fname}.hotspot.nooverlap.vcf.gz
  tabix ${gatk}
  tabix ${sam}
  tabix ${ss}
  java -Xmx32g -jar \$GATK_JAR -R ${index_path}/hisat_genome.fa -T CombineVariants --filteredrecordsmergetype KEEP_UNCONDITIONAL --variant:gatk ${gatk} --variant:sam ${sam} --variant:ssvar ${ss} --variant:hotspot ${fname}.hotspot.nooverlap.vcf.gz -genotypeMergeOptions PRIORITIZE -priority sam,ssvar,gatk,hotspot -o ${fname}.int.vcf
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
  bedtools multiinter -i ${gatk} ${sam} ${ss} -names gatk sam ssvar platypus |cut -f 1,2,3,5 | bedtools sort -i stdin | bedtools merge -c 4 -o distinct >  ${fname}_integrate.bed
  java -Xmx32g -jar \$GATK_JAR -R ${index_path}/hisat_genome.fa -T CombineVariants --filteredrecordsmergetype KEEP_UNCONDITIONAL -genotypeMergeOptions PRIORITIZE --variant:gatk ${gatk} --variant:sam ${sam} --variant:ssvar ${ss} -priority sam,ssvar,gatk -o ${fname}.int.vcf
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
  file("${fname}.stats.txt") into stats
  file("${fname}.statplot*") into plotstats
  when:
  params.variant == "detect"
  script:
  if (params.genome == '/project/shared/bicf_workflow_ref/GRCh38')
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
  """
  else
  """
  module load snpeff/4.2
  java -Xmx10g -jar \$SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c \$SNPEFF_HOME/snpEff.config ${snpeff_vers} ${unionvcf} |bgzip > ${fname}.annot.vcf.gz
  tabix ${fname}.annot.vcf.gz
  bcftools stats ${fname}.annot.vcf.gz > ${fname}.stats.txt
  plot-vcfstats -s -p ${fname}.statplot ${fname}.stats.txt
  """
}
//
// Read summarization with subread
//
process featurect {

  publishDir "$params.output", mode: 'copy'
  cpus 32

  input:
  set pair_id, file(dbam) from deduped1
  file gtf_file
  output:
  file("${pair_id}.cts")  into counts
  """
  module load module subread/1.5.0-intel
  featureCounts -s params.stranded -T 30 -p -g gene_name -a ${gtf_file} -o ${pair_id}.cts ${dbam}
  """
}

//
// Assemble transcripts with stringtie
//

process stringtie {

  publishDir "$params.output", mode: 'copy'
  cpus 32

  input:
  set pair_id, file(dbam) from deduped2
  file stringtie_gtf
  output:
  file("${pair_id}_stringtie") into strcts
  file("${pair_id}.fpkm.txt") into fpkm
  """
  module load stringtie/1.1.2-intel
  mkdir ${pair_id}_stringtie
  cd ${pair_id}_stringtie
  stringtie ../${dbam} -p 30 -G ../${stringtie_gtf} -B -e -o denovo.gtf -A ../${pair_id}.fpkm.txt
  """
}

process statanal {
   errorStrategy 'ignore'
   publishDir "$params.output", mode: 'copy'
   input:
   file count_file from counts.toList()
   file design_file name 'design.txt'
   file genenames
   file geneset name 'geneset.gmt'
   output:
   file "*" into txtfiles
 
   script:
   if (params.dea == 'detect')
   """
   module load R/3.2.1-intel
   perl $baseDir/scripts/concat_cts.pl -o ./ *.cts
   cp design.txt design.shiny.txt
   cp geneset.gmt geneset.shiny.gmt
   Rscript  $baseDir/scripts/dea.R
   perl $baseDir/scripts/concat_edgeR.pl *.edgeR.txt
   """
   else
   """
   perl $baseDir/scripts/concat_cts.pl -o ./ *.cts
   """
}

process buildrda {
   publishDir "$params.output", mode: 'copy'
   input:
   file stringtie_dir from strcts.toList()
   file design_file name 'design.txt'
   output:
   file "bg.rda" into rdafiles
   when:
   params.dea == 'detect'
   script:
   """
   module load R/3.2.1-intel
   Rscript $baseDir/scripts/build_ballgown.R *_stringtie
   """
 }
