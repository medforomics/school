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

process svcall {
  publishDir "$params.output", mode: 'copy'
  errorStrategy 'ignore'
  input:
  set tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from dellybam

  output:
  file("${tid}_${nid}.delly.vcf.gz") into dellyvcf
  file("${tid}_${nid}.novobreak.vcf.gz") into novovcf
  file("${tid}_${nid}.sssv.sv.vcf.gz") into lumpyvcf
  file("${tid}_${nid}.sv.vcf.gz") into svintvcf
  file("${tid}_${nid}.sv.annot.txt") into svannot
  script:
  """
  source /etc/profile.d/modules.sh
  module load novoBreak/v1.1.3 delly2/v0.7.7-multi bcftools/intel/1.3 samtools/intel/1.3 bedtools/2.25.0 speedseq/20160506 snpeff/4.2 vcftools/0.1.14
  mkdir temp
  perl $baseDir/scripts/make_delly_sample.pl ${tid} ${nid}
  delly2 call -t BND -o delly_translocations.bcf -q 30 -g ${reffa} ${tumor} ${normal}
  delly2 call -t DUP -o delly_duplications.bcf -q 30 -g ${reffa} ${tumor} ${normal}
  delly2 call -t INV -o delly_inversions.bcf -q 30 -g ${reffa} ${tumor} ${normal}
  delly2 call -t DEL -o delly_deletion.bcf -q 30 -g ${reffa} ${tumor} ${normal}
  delly2 call -t INS -o delly_insertion.bcf -q 30 -g ${reffa} ${tumor} ${normal}
  delly2 filter -t BND -o  delly_tra.bcf -f somatic -s samples.tsv delly_translocations.bcf
  delly2 filter -t DUP -o  delly_dup.bcf -f somatic -s samples.tsv delly_translocations.bcf
  delly2 filter -t INV -o  delly_inv.bcf -f somatic -s samples.tsv delly_translocations.bcf
  delly2 filter -t DEL -o  delly_del.bcf -f somatic -s samples.tsv delly_translocations.bcf
  delly2 filter -t INS -o  delly_ins.bcf -f somatic -s samples.tsv delly_translocations.bcf
  bcftools concat -a -O v delly_dup.bcf delly_inv.bcf delly_tra.bcf delly_del.bcf delly_ins.bcf| vcf-sort -t temp > ${tid}_${nid}.delly.vcf
  perl $baseDir/scripts/vcf2bed.sv.pl ${tid}_${nid}.delly.vcf |sort -T temp -V -u -k 1,1 -k 2,2n > delly.bed
  bgzip ${tid}_${nid}.delly.vcf
  tabix ${tid}_${nid}.delly.vcf.gz
  bcftools view -O z -o delly.vcf.gz -s ${tid} ${tid}_${nid}.delly.vcf.gz
  run_novoBreak.sh /cm/shared/apps/novoBreak/novoBreak_distribution_v1.1.3rc ${reffa} ${tumor} ${normal} \$SLURM_CPUS_ON_NODE
  perl $baseDir/scripts/vcf2bed.sv.pl novoBreak.pass.flt.vcf |sort -T temp -V -k 1,1 -k 2,2n > novobreak.bed
  mv novoBreak.pass.flt.vcf ${tid}_${nid}.novobreak.vcf
  bgzip ${tid}_${nid}.novobreak.vcf
  sambamba sort --tmpdir=./ -t \$SLURM_CPUS_ON_NODE -n -o tumor.namesort.bam ${tumor}
  sambamba sort --tmpdir=./ -t \$SLURM_CPUS_ON_NODE -n -o normal.namesort.bam ${normal}
  sambamba view -h tumor.namesort.bam | samblaster -M -a --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 -d discordants.sam -s splitters.sam > temp.sam
  gawk '{ if (\$0~"^@") { print; next } else { \$10="*"; \$11="*"; print } }' OFS="\\t" splitters.sam | samtools  view -S -b - | samtools sort -o tumor.splitters.bam -
  gawk '{ if (\$0~"^@") { print; next } else { \$10="*"; \$11="*"; print } }' OFS="\\t" discordants.sam | samtools  view -S  -b - | samtools sort -o tumor.discordants.bam -
  sambamba view -h normal.namesort.bam | samblaster -M -a --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 -d discordants.sam -s splitters.sam > temp.sam
  gawk '{ if (\$0~"^@") { print; next } else { \$10="*"; \$11="*"; print } }' OFS="\\t" splitters.sam | samtools  view -S -b - | samtools sort -o normal.splitters.bam -
  gawk '{ if (\$0~"^@") { print; next } else { \$10="*"; \$11="*"; print } }' OFS="\\t" discordants.sam | samtools  view -S  -b - | samtools sort -o normal.discordants.bam -
  speedseq sv -t \$SLURM_CPUS_ON_NODE -o ${tid}_${nid}.sssv -R ${reffa} -B ${normal},${tumor} -D normal.discordants.bam,tumor.discordants.bam -S normal.splitters.bam,tumor.splitters.bam -x ${index_path}/exclude_alt.bed
  java -jar \$SNPEFF_HOME/SnpSift.jar filter "GEN[0].SU < 1 & GEN[1].SU > 2" ${tid}_${nid}.sssv.sv.vcf.gz > lumpy.vcf
  perl $baseDir/scripts/vcf2bed.sv.pl lumpy.vcf |sort -T temp -V -u -k 1,1 -k 2,2n > lumpy.bed
  bcftools view -O z -o sssv.vcf.gz -s ${tid} ${tid}_${nid}.sssv.sv.vcf.gz 
  bedtools multiinter -cluster -header -names novobreak delly lumpy -i novobreak.bed delly.bed lumpy.bed > sv.intersect.bed 
  grep novobreak sv.intersect.bed |cut -f 1,2,3 |sort -V -k 1,1 -k 2,2n |grep -v start | bedtools intersect -header -b stdin -a ${tid}_${nid}.novobreak.vcf.gz  | perl -p -e 's/SPIKEIN/${tid}/' |bgzip > t1.vcf.gz
  grep delly sv.intersect.bed |cut -f 1,2,3 |sort -V -k 1,1 -k 2,2n |grep -v 'start' |grep -v 'novobreak' | bedtools intersect -header -b stdin -a delly.vcf.gz |bgzip > t2.vcf.gz
  grep lumpy sv.intersect.bed |cut -f 1,2,3 |sort -V -k 1,1 -k 2,2n |grep -v 'start' |grep -v 'delly' |grep -v 'novobreak' | bedtools intersect -header -b stdin -a sssv.vcf.gz |bgzip > t3.vcf.gz
  vcf-concat t1.vcf.gz t2.vcf.gz t3.vcf.gz |vcf-sort -t temp > ${tid}_${nid}.sv.vcf
  perl $baseDir/scripts/vcf2bed.sv.pl ${tid}_${nid}.sv.vcf |sort -V -k 1,1 -k 2,2n | grep -v 'alt' |grep -v 'random' |uniq > svs.bed
  bedtools intersect -header -wb -a svs.bed -b ${index_path}/gencode.exons.bed > exonoverlap_sv.txt
  bedtools intersect -v -header -wb -a svs.bed -b ${index_path}/gencode.exons.bed | bedtools intersect -header -wb -a stdin -b ${index_path}/gencode.genes.chr.bed > geneoverlap_sv.txt
  perl $baseDir/scripts/annot_sv.pl -r ${index_path} -i ${tid}_${nid}.sv.vcf
  bgzip ${tid}_${nid}.sv.vcf
  """
}
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
  cut -f 1 ${index_path}/genomefile.5M.txt | parallel --delay 2 -j 10 "java -Xmx20g -jar \$GATK_JAR -R ${reffa} -D ${dbsnp} -T MuTect2 -stand_call_conf 30 -stand_emit_conf 10.0 -A FisherStrand -A QualByDepth -A VariantType -A DepthPerAlleleBySample -A HaplotypeScore -A AlleleBalance -I:tumor ${tumor} -I:normal ${normal} --cosmic ${cosmic} -o ${tid}.{}.mutect.vcf -L {}"
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

ssvcf .phase(mutectvcf)
      .map {p,q -> [p[0],p[1],q[1]]}
      .set { twovcf }
twovcf .phase(varscanvcf)
      .map {p,q -> [p[0],p[1],p[2],q[1]]}
      .set { threevcf }
threevcf .phase(virmidvcf)
      .map {p,q -> [p[0],p[1],p[2],p[3],q[1]]}
      .set { fourvcf }
fourvcf .phase(shimmervcf)
	.map {p,q -> [p[0],p[1],p[2],p[3],p[4],q[1]]}
      	.set { vcflist }

process integrate {
  errorStrategy 'ignore'
  //publishDir "$params.output", mode: 'copy'

  input:
  set fname,file(ss),file(mutect),file(vscan),file(virmid),file(shimmer) from vcflist
  
  output:
  set fname,file("${fname}.union.vcf.gz") into union
  script:
  """
  source /etc/profile.d/modules.sh
  module load gatk/3.5 python/2.7.x-anaconda bedtools/2.25.0 snpeff/4.2 bcftools/intel/1.3 samtools/intel/1.3
  module load vcftools/0.1.14
  perl $baseDir/scripts/somatic_unionize_vcf.pl -r ${index_path} ${ss} ${mutect} ${shimmer} ${vscan} ${virmid}
  sh integrate.sh
  perl $baseDir/scripts/somatic_uniform_integrated_vcf.pl ${fname}.temp.vcf
  bgzip ${fname}.union.vcf
  """
}

union .phase(pairnames)
      .map {p,q -> [p[0],p[1],q[1],q[2]]}
      .set { union2 }
      
process annot {
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'

  input:
  set fname,unionvcf,tid,nid from union2
  
  output:
  file("${fname}.annot.vcf.gz") into annot

  script:
  """
  source /etc/profile.d/modules.sh
  module load python/2.7.x-anaconda bedtools/2.25.0 snpeff/4.2 bcftools/intel/1.3 samtools/intel/1.3 git/v2.5.3
  tabix ${unionvcf}
  bcftools annotate -Oz -a ${index_path}/ExAC.vcf.gz -o ${fname}.exac.vcf.gz --columns CHROM,POS,AC_Het,AC_Hom,AC_Hemi,AC_Adj,AN_Adj,AC_POPMAX,AN_POPMAX,POPMAX ${unionvcf}
  tabix ${fname}.exac.vcf.gz 
  java -Xmx10g -jar \$SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c \$SNPEFF_HOME/snpEff.config ${snpeff_vers} ${fname}.exac.vcf.gz | java -jar \$SNPEFF_HOME/SnpSift.jar annotate -id ${index_path}/dbSnp.vcf.gz -  | java -jar \$SNPEFF_HOME/SnpSift.jar annotate -info CLNSIG,CLNDSDB,CLNDSDBID,CLNDBN,CLNREVSTAT,CLNACC ${index_path}/clinvar.vcf.gz - | java -jar \$SNPEFF_HOME/SnpSift.jar annotate -info CNT ${index_path}/cosmic.vcf.gz - | java -Xmx10g -jar \$SNPEFF_HOME/SnpSift.jar dbnsfp -v -db ${index_path}/dbNSFP.txt.gz - | bgzip > ${fname}.annot.vcf.gz
  tabix ${fname}.annot.vcf.gz
  """
}
