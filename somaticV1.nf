#!/usr/bin/env nextflow

// Default parameter values to run tests
params.input = './'
params.output = './'

params.bams="$params.input/*.bam"
params.design="$params.input/design.txt"

params.genome="/project/shared/bicf_workflow_ref/GRCh38"
params.targetpanel="$params.genome/UTSWV2.bed"

dbsnp="$params.genome/dbSnp.vcf.gz"
cosmic="$params.genome/cosmic.vcf.gz"

design_file = file(params.design)
bams=file(params.bams)

reffa=file("$params.genome/genome.fa")
index_path = file(params.genome)
ncmconf = file("$params.genome/ncm.conf")
target_panel = file(params.targetpanel)
dbsnp=file(dbsnp)

strelkaconfig="/cm/shared/apps/strelka/1.0.15/etc/strelka_config_bwa_default.ini"
confstrelka=file(strelkaconfig)

snpeff_vers = 'GRCh38.82';
if (params.genome == '/project/shared/bicf_workflow_ref/GRCm38') {
   snpeff_vers = 'GRCm38.82';
}
if (params.genome == '/project/shared/bicf_workflow_ref/GRCh37') {
   snpeff_vers = 'GRCh37.75';
}

def fileMap = [:]

bams.each {
    final fileName = it.getFileName().toString()
    prefix = fileName.lastIndexOf('/')
    fileMap[fileName] = it
}
def oribam = []
def tarbam = []
new File(params.design).withReader { reader ->
    def hline = reader.readLine()
    def header = hline.split("\t")
    tidx = header.findIndexOf{it == 'TumorID'};
    nidx = header.findIndexOf{it == 'NormalID'};
    oneidx = header.findIndexOf{it == 'TumorBAM'};
    twoidx = header.findIndexOf{it == 'NormalBAM'};
    totidx = header.findIndexOf{it == 'TumorOntargetBAM'};
    notidx = header.findIndexOf{it == 'NormalOntargetBAM'};

    if (twoidx == -1) {
       twoidx = oneidx
       }      
    while (line = reader.readLine()) {
    	   def row = line.split("\t")
	   if (fileMap.get(row[oneidx]) != null) {
	      oribam << tuple(row[tidx],row[nidx],fileMap.get(row[oneidx]),fileMap.get(row[twoidx]))
	      tarbam << tuple(row[tidx],row[nidx],fileMap.get(row[totidx]),fileMap.get(row[notidx]))
	   }
} 
}

if( ! oribam) { error "Didn't match any input files with entries in the design file" }
if( ! tarbam) { error "Didn't match any input files with entries in the design file" }

process indexoribams {
  errorStrategy 'ignore'
  input:
  set tid,nid,file(tumor),file(normal) from oribam
  output:
  set tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into dellybam
  set tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into checkbams

  script:
  """
  source /etc/profile.d/modules.sh
  module load speedseq/20160506 samtools/intel/1.3
  sambamba index -t \$SLURM_CPUS_ON_NODE ${tumor}
  sambamba index -t \$SLURM_CPUS_ON_NODE ${normal}
  """
}
process indextarbams {
  errorStrategy 'ignore'
  input:
  set tid,nid,file(tumor),file(normal) from tarbam
  output:
  set tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into mutectbam
  set tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into ssbam
  set tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into shimmerbam
  set tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into vscanbam
  set tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into virmidbam
  set val("${tid}_${nid}"),tid,nid into pairnames
  script:
  """
  source /etc/profile.d/modules.sh
  module load speedseq/20160506 samtools/intel/1.3 gatk/3.5
  sambamba index -t \$SLURM_CPUS_ON_NODE ${tumor}
  sambamba index -t \$SLURM_CPUS_ON_NODE ${normal}
  """
}
process checkmates {
  publishDir "$params.output", mode: 'copy'
  errorStrategy 'ignore'
  input:
  set tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from checkbams
  file(conf) from ncmconf
  output:
  file("${tid}_${nid}*") into checkmateout
  script:
  """
  source /etc/profile.d/modules.sh
  module load python/2.7.x-anaconda bcftools/intel/1.3 samtools/intel/1.3 git/v2.5.3
  python /project/shared/bicf_workflow_ref/seqprg/NGSCheckMate/ncm.py -B -d ./ -bed ${index_path}/NGSCheckMate.bed -O ./ -N ${tid}_${nid}
  perl $baseDir/scripts/sequenceqc_somatic.pl -r ${index_path} -i ${tid}_${nid}_all.txt -o ${tid}_${nid}.sequence.stats.txt
  """
}

// process svcall {
//   publishDir "$params.output", mode: 'copy'
//   errorStrategy 'ignore'
//   input:
//   set tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from dellybam

//   output:
//   file("${tid}_${nid}.delly.vcf.gz") into dellyvcf
//   file("${tid}_${nid}.sssv.sv.vcf.gz") into lumpyvcf
//   file("${tid}_${nid}.sv.annot.vcf.gz") into svintvcf
//   file("${tid}_${nid}.sv.annot.txt") into svannot
//   file("${tid}_${nid}.sv.annot.genefusion.txt") into gfannot
//   when:
//   params.callsvs == "detect"
//   script:
//   """
//   source /etc/profile.d/modules.sh
//   module load novoBreak/v1.1.3 delly2/v0.7.7-multi bcftools/intel/1.3 samtools/intel/1.3 bedtools/2.25.0 speedseq/20160506 snpeff/4.2 vcftools/0.1.14
//   perl $baseDir/scripts/make_delly_sample.pl ${tid} ${nid}
//   bash $baseDir/process_scripts/variants/svcalling.sh -r ${index_path} -p ${tid}_${nid} -b ${tumor} -n ${normal} -k ${tid}
//   """
// }

process sstumor {
  errorStrategy 'ignore'
  //publishDir "$params.output", mode: 'copy'
  input:
  set tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from ssbam
  output:
  set val("${tid}_${nid}"), file("${tid}_${nid}.sspanel.vcf.gz") into ssvcf
  script:
  """
  source /etc/profile.d/modules.sh
  module load python/2.7.x-anaconda bedtools/2.25.0 snpeff/4.2 speedseq/20160506 bcftools/intel/1.3 vcftools/0.1.14
  speedseq somatic -q 10 -w ${target_panel} -t \$SLURM_CPUS_ON_NODE -o ${tid}.sssom ${reffa} ${normal} ${tumor}
  vcf-annotate -H -n --fill-type ${tid}.sssom.vcf.gz | java -jar \$SNPEFF_HOME/SnpSift.jar filter --pass '((QUAL >= 10) & (GEN[*].DP >= 10))' | perl -pe 's/TUMOR/${tid}/' | perl -pe 's/NORMAL/${nid}/g' |bgzip > ${tid}_${nid}.sspanel.vcf.gz
  """
}
process mutect {
  errorStrategy 'ignore'
  //publishDir "$params.output", mode: 'copy'

  input:
  set tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from mutectbam

  output:
  set val("${tid}_${nid}"),file("${tid}_${nid}.pmutect.vcf.gz") into mutectvcf
  script:
  """
  source /etc/profile.d/modules.sh
  module load parallel python/2.7.x-anaconda gatk/3.5  bcftools/intel/1.3 bedtools/2.25.0 snpeff/4.2 vcftools/0.1.14
  cut -f 1 ${index_path}/genomefile.5M.txt | parallel --delay 2 -j 6 "java -Xmx20g -jar \$GATK_JAR -R ${reffa} -D ${dbsnp} -T MuTect2 -stand_call_conf 30 -stand_emit_conf 10.0 -A FisherStrand -A QualByDepth -A VariantType -A DepthPerAlleleBySample -A HaplotypeScore -A AlleleBalance -I:tumor ${tumor} -I:normal ${normal} --cosmic ${cosmic} -o ${tid}.{}.mutect.vcf -L {}"
  vcf-concat ${tid}*.vcf | vcf-sort | vcf-annotate -n --fill-type | java -jar \$SNPEFF_HOME/SnpSift.jar filter -p '((FS <= 60) & GEN[*].DP >= 10)' | perl -pe 's/TUMOR/${tid}/' | perl -pe 's/NORMAL/${nid}/g' |bgzip > ${tid}_${nid}.pmutect.vcf.gz
  """
}
process varscan {
  errorStrategy 'ignore'
  //publishDir "$params.output", mode: 'copy'
  input:
  set tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from vscanbam
  output:
  set val("${tid}_${nid}"),file("${tid}_${nid}.varscan.vcf.gz") into varscanvcf
  script:
  """
  source /etc/profile.d/modules.sh
  module load python/2.7.x-anaconda bedtools/2.25.0 snpeff/4.2 bcftools/intel/1.3 samtools/intel/1.3 VarScan/2.4.2 speedseq/20160506 vcftools/0.1.14
  sambamba mpileup --tmpdir=./ -L ${target_panel} -t \$SLURM_CPUS_ON_NODE ${tumor} --samtools "-C 50 -f ${reffa}"  > t.mpileup
  sambamba mpileup --tmpdir=./ -L ${target_panel} -t \$SLURM_CPUS_ON_NODE ${normal} --samtools "-C 50 -f ${reffa}"  > n.mpileup
  VarScan somatic n.mpileup t.mpileup ${tid}.vscan --output-vcf 1
  VarScan copynumber n.mpileup t.mpileup ${tid}.vscancnv 
  vcf-concat ${tid}.vscan*.vcf | vcf-sort | vcf-annotate -n --fill-type -n | java -jar \$SNPEFF_HOME/SnpSift.jar filter '((exists SOMATIC) & (GEN[*].DP >= 10))' | perl -pe 's/TUMOR/${tid}/' | perl -pe 's/NORMAL/${nid}/g' |bedtools intersect -header -a stdin -b ${target_panel} |bgzip > ${tid}_${nid}.varscan.vcf.gz
  """
}
process shimmer {
  errorStrategy 'ignore'
  //publishDir "$params.output", mode: 'copy'
  input:
  set tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from shimmerbam
  output:
  set val("${tid}_${nid}"), file("${tid}_${nid}.shimmer.vcf.gz") into shimmervcf
  script:
  """
  source /etc/profile.d/modules.sh
  module load python/2.7.x-anaconda bedtools/2.25.0 snpeff/4.2 bcftools/intel/1.3  shimmer/0.1.1 vcftools/0.1.14
  shimmer.pl --minqual 25 --ref ${reffa} ${normal} ${tumor} --outdir shimmer 2> shimmer.err
  perl $baseDir/scripts/add_readct_shimmer.pl
  vcf-annotate -n --fill-type shimmer/somatic_diffs.readct.vcf | java -jar \$SNPEFF_HOME/SnpSift.jar filter '(GEN[*].DP >= 10)' | perl -pe 's/TUMOR/${tid}/' | perl -pe 's/NORMAL/${nid}/g' | bedtools intersect -header -a stdin -b ${target_panel} | bgzip > ${tid}_${nid}.shimmer.vcf.gz
  """
}

process virmid {
  errorStrategy 'ignore'
  //publishDir "$params.output", mode: 'copy'
  input:
  set tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from virmidbam
  output:
  set val("${tid}_${nid}"), file("${tid}_${nid}.virmid.vcf.gz") into virmidvcf
  script:
  """
  source /etc/profile.d/modules.sh
  module load python/2.7.x-anaconda bedtools/2.25.0 snpeff/4.2 virmid/1.2 vcftools/0.1.14
  virmid -R ${reffa} -D ${tumor} -N ${normal} -s $cosmic -t \$SLURM_CPUS_ON_NODE -M 2000 -c1 10 -c2 10
  perl $baseDir/scripts/addgt_virmid.pl ${tumor}.virmid.som.passed.vcf
  perl $baseDir/scripts/addgt_virmid.pl ${tumor}.virmid.loh.passed.vcf
  vcf-concat *gt.vcf | vcf-sort | vcf-annotate -n --fill-type -n | java -jar \$SNPEFF_HOME/SnpSift.jar filter '((NDP >= 10) & (DDP >= 10))' | perl -pe 's/TUMOR/${tid}/' | perl -pe 's/NORMAL/${nid}/g' |bedtools intersect -header -a stdin -b ${target_panel} |bgzip > ${tid}_${nid}.virmid.vcf.gz
  """
}

Channel
  .empty()
  .mix(ssvcf,mutectvcf,varscanvcf,virmidvcf,shimmervcf)
  .groupTuple(by:0)
  .into { vcflist}


process integrate {
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'
  input:
  set subjid,file(vcf) from vcflist
  output:
  set subjid,file("${subjid}.union.vcf.gz") into union
  file("${subjid}.annot.vcf.gz") into annotvcf
  script:
  """
  bash $baseDir/process_scripts/variants/union.sh -r $index_path -p $subjid
  bash $baseDir/process_scripts/variants/annotvcf.sh -p $subjid -r $index_path -v ${subjid}.union.vcf.gz
  """
}
