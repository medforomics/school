#!/usr/bin/env nextflow
	
// Default parameter values to run tests
params.input = './fastq'
params.output = './'

params.fastqs="$params.input/*.fastq.gz"
params.design="$params.input/design.txt"

params.genome="/project/shared/bicf_workflow_ref/human/grch38_cloud/rnaref/"
params.markdups="skip"
params.stranded="0"
params.pairs="pe"
params.align = 'hisat'
params.fusion = 'detect'
params.dea = 'skip'
params.cancer="detect"
params.variant = "detect"
params.bamct = "detect"
params.gatk='skip'

params.version='v4'

design_file = file(params.design)
fastqs=file(params.fastqs)
design_file = file(params.design)
gtf_file = file("$params.genome/gencode.gtf")
genenames = file("$params.genome/genenames.txt")
dbsnp="$params.genome/dbSnp.vcf.gz"
indel="$params.genome/GoldIndels.vcf.gz"
knownindel=file(indel)
dbsnp=file(dbsnp)

snpeff_vers = 'GRCh38.86';

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
process rtrim {
  executor 'local'
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
  bash $baseDir/process_scripts/alignment/starfusion.sh -p ${pair_id} -r ${index_path} -a ${fq1} -b ${fq2} -m trinity -f
  """
}
process ralign {
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid/$pair_id", mode: 'copy'

  input:
  set subjid,pair_id, file(f1), file(f2) from trimread
  output:
  set subjid,pair_id, file("${pair_id}.bam") into aligned
  set subjid,pair_id, file("${pair_id}.bam") into ctbams
  set subjid,pair_id, file("${pair_id}.bam"),file("${pair_id}.alignerout.txt") into aligned2
  set subjid,pair_id, file("${pair_id}.bam") into deduped1
  set subjid,pair_id, file("${pair_id}.bam"),file("${pair_id}.bam.bai") into fbbam

  script:
  """
  bash $baseDir/process_scripts/alignment/rnaseqalign.sh -a $params.align -p $pair_id -r $index_path -x $f1 -y $f2 $alignopts
  """
}
process bamct {
  executor 'local'
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
  /project/shared/bicf_workflow_ref/seqprg/bam-readcount/bin/bam-readcount -w 0 -q 0 -b 25 -f ${index_path}/genome.fa $rbam > ${pair_id}.bamreadct.txt
  """
}

process alignqc {
  executor 'local'
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
  perl $baseDir/scripts/sequenceqc_rna.pl -r ${index_path} -e $params.version *.flagstat.txt
  """
}

// Identify duplicate reads with Picard


process geneabund {
  errorStrategy 'ignore'
  executor 'local'
  publishDir "$params.output/$subjid/$pair_id", mode: 'copy'
  input:
  set subjid,pair_id, file(sbam) from deduped1
  output:
  file("${pair_id}.cts")  into counts
  file("${pair_id}_stringtie") into strcts
  file("${pair_id}.fpkm.txt") into fpkm
  """
  source /etc/profile.d/modules.sh
  bash $baseDir/process_scripts/diff_exp/geneabundance.sh -s $params.stranded -g ${gtf_file} -p ${pair_id} -b ${sbam} -f 1
  """
}

process fb {
  queue '32GB'
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid", mode: 'copy'

  input:
  set subjid,$pair_id,file(gbam),file(gidx) from fbbam
  output:
  set subjid,file("${subjid}.rna.vcf.gz") into fbvcf
  script:
  """
  bash $baseDir/process_scripts/variants/germline_vc.sh -r $index_path -p $subjid -a fb
  bash $baseDir/process_scripts/variants/uni_norm_annot.sh -g $snpeff_vers -r $index_path -p ${subjid}.fb -v ${subjid}.fb.vcf.gz 
  mv ${subjid}.fb.vcf.gz ${subjid}_${params.projectid}.rna.vcf.gz
  """
}
