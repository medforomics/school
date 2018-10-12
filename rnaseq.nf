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
params.bamct = "detect"


design_file = file(params.design)
fastqs=file(params.fastqs)
design_file = file(params.design)
gtf_file = file("$params.genome/gencode.gtf")
stringtie_gtf=file("$params.genome/gencode.hisat.gtf")
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
alignopts = ''
if (params.markdups == 'fgbio_umi') {
   alignopts = '-u'
}
// params genome is the directory
// base name for the index is always genome
index_path = file(params.genome)
index_name = "genome"

// Pair handling, helper function taken from rnatoy
// which is covered by the GNU General Public License v3
// https://github.com/nextflow-io/rnatoy/blob/master/main.nf

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
    tidx = header.findIndexOf{it == 'SampleID'};
    oneidx = header.findIndexOf{it == 'FqR1'};
    twoidx = header.findIndexOf{it == 'FqR2'};
    fidx = header.findIndexOf{it == 'FamilyID'};
    sidx = header.findIndexOf{it == 'SubjectID'};
    if (sidx == -1) {
       sidx = tidx
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
	      read << tuple(row[fidx],row[tidx],fileMap.get(row[oneidx]),fileMap.get(row[twoidx]))
	   }
	  
} 
}
if( ! read) { error "Didn't match any input files with entries in the design file" }

//
// Trim raw reads using trimgalore
//
process trim {
  errorStrategy 'ignore'
  input:
  set subjid, pair_id, file(read1), file(read2) from read
  output:
  set subjid, pair_id, file("${pair_id}.trim.R1.fastq.gz"),file("${pair_id}.trim.R2.fastq.gz") into trimread
  set subjid, pair_id, file("${pair_id}.trim.R1.fastq.gz"),file("${pair_id}.trim.R2.fastq.gz") into fusionfq
  script:
  """
  bash $baseDir/process_scripts/preproc_fastq/trimgalore.sh -p ${pair_id} -a ${read1} -b ${read2}
  """
}
process starfusion {
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid/$pair_id", mode: 'copy'
  input:
  set subjid,pair_id, file(fq1), file(fq2) from fusionfq
  output:
  file("${pair_id}*txt") into fusionout
  when:
  params.fusion == 'detect' && params.pairs == 'pe'
  script:
  """
  bash $baseDir/process_scripts/alignment/starfusion.sh -p ${pair_id} -r ${index_path} -a ${fq1} -b ${fq2} -m trinity
  """
}
process align {
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid/$pair_id", mode: 'copy'

  input:
  set subjid,pair_id, file(f1), file(f2) from trimread
  output:
  set subjid,pair_id, file("${pair_id}.bam") into aligned
  set subjid,pair_id, file("${pair_id}.bam") into ctbams
  set subjid,pair_id, file("${pair_id}.bam"),file("${pair_id}.alignerout.txt") into aligned2

  script:
  """
  bash $baseDir/process_scripts/alignment/rnaseqalign.sh -a $params.align -p $pair_id -r $index_path -x $f1 -y $f2 $alignopts
  """
}
process bamct {
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid/$pair_id", mode: 'copy'
  input:
  set subjid,pair_id, file(rbam) from ctbams
  output:
  file("${pair_id}.bamreadct.txt") into ctreads
  when:
  params.bamct == "detect"
  script:
  """
  source /etc/profile.d/modules.sh
  module load samtools/1.6
  samtools index $rbam
  ${index_path}/../seqprg/bam-readcount/bin/bam-readcount -w 0 -q 0 -b 25 -f ${index_path}/hisat_genome.fa $rbam > ${pair_id}.bamreadct.txt
  """
}

process alignqc {
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid/$pair_id", mode: 'copy'

  input:
  set subjid,pair_id, file(bam), file(hsout) from aligned2
  
  output:
  set file("${pair_id}_fastqc.zip"),file("${pair_id}_fastqc.html") into fastqc
  file("${pair_id}.sequence.stats.txt") into alignstats
  script:
  """
  source /etc/profile.d/modules.sh
  module load git/gcc/v2.12.2
  bash $baseDir/process_scripts/alignment/bamqc.sh -p ${pair_id} -b ${bam} -n rna
  perl $baseDir/scripts/sequenceqc_rnaseq.pl -r ${index_path} *.flagstat.txt
  """
}

// Identify duplicate reads with Picard

process markdups {
  //publishDir "$params.output/$subjid/$pair_id", mode: 'copy'

  input:
  set subjid,pair_id, file(sbam) from aligned
  output:
  set subjid,pair_id, file("${pair_id}.dedup.bam") into deduped1
  set subjid,pair_id, file("${pair_id}.dedup.bam") into deduped2
  script:
  """
  bash $baseDir/process_scripts/alignment/markdups.sh -a $params.markdups -b $sbam -p $pair_id
  """
}

process geneabund {
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid/$pair_id", mode: 'copy'
  input:
  set subjid,pair_id, file(sbam) from deduped1
  file gtf_file
  output:
  file("${pair_id}.cts")  into counts
  file("${pair_id}_stringtie") into strcts
  file("${pair_id}.fpkm.txt") into fpkm
  file('*cbioportal.txt')
  """
  source /etc/profile.d/modules.sh
  bash $baseDir/process_scripts/diff_exp/geneabundance.sh -s $params.stranded -g ${gtf_file} -p ${pair_id} -b ${sbam}
  perl $baseDir/process_scripts/genect_rnaseq/cBioPortal_documents.pl -p $pair_id -l ${pair_id}.cts -f ${pair_id}.fpkm.txt
  """
}
process gatkbam {
  errorStrategy 'ignore'
  publishDir "$params.output/", mode: 'copy'

  input:
  set subjid,pair_id, file(rbam) from deduped2

  output:
  set file("${pair_id}.final.bam"),file("${pair_id}.final.bai") into gatkbam

  script:
  """
  bash $baseDir/process_scripts/variants/gatkrunner.sh -a gatkbam_rna -b $rbam -r ${index_path}/hisat_index -p $pair_id
  """
}
