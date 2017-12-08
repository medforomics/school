#!/usr/bin/env nextflow

params.input = './analysis'
params.output = './analysis'

params.bams="$params.input/*.bam"
params.design="$params.input/design.txt"

params.genome="/project/shared/bicf_workflow_ref/GRCh38"
params.targetpanel="$params.genome/UTSWV2.bed"
params.cancer="detect"

dbsnp_file="$params.genome/dbSnp.vcf.gz"
indel="$params.genome/GoldIndels.vcf.gz"
cosmic="$params.genome/cosmic.vcf.gz"
reffa=file("$params.genome/genome.fa")

design_file = file(params.design)
bams=file(params.bams)

dbsnp=file(dbsnp_file)
knownindel=file(indel)
index_path = file(params.genome)
capture_bed = file(params.targetpanel)

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
    tidx = header.findIndexOf{it == 'SampleID'};
    oneidx = header.findIndexOf{it == 'BAM'};
    taridx = header.findIndexOf{it == 'OntargetBAM'};
    while (line = reader.readLine()) {
    	   def row = line.split("\t")
	   if (fileMap.get(row[oneidx]) != null) {
	      oribam << tuple(row[tidx],fileMap.get(row[oneidx]))
	      tarbam << tuple(row[tidx],fileMap.get(row[taridx]))
	   }
	  
} 
}

if( ! oribam) { error "Didn't match any input files with entries in the design file" }
if( ! tarbam) { error "Didn't match any input files with entries in the design file" }

process indexoribams {
  errorStrategy 'ignore'
  input:
  set tid,file(tumor) from oribam
  output:
  set tid,file(tumor),file("${tumor}.bai") into ssbam
  set tid,file(tumor),file("${tumor}.bai") into svbam
  script:
  """
  source /etc/profile.d/modules.sh
  module load speedseq/20160506 samtools/intel/1.3
  sambamba index -t \$SLURM_CPUS_ON_NODE ${tumor}
  """
}
process indexbams {
  input:
  set tid,file(tumor) from tarbam
  output:
  set tid,file(tumor),file("${tumor}.bai") into gatkbam
  set tid,file(tumor),file("${tumor}.bai") into sambam
  set tid,file(tumor),file("${tumor}.bai") into hsbam
  set tid,file(tumor),file("${tumor}.bai") into platbam

  script:
  """
  source /etc/profile.d/modules.sh
  module load speedseq/20160506 samtools/intel/1.3
  sambamba index -t \$SLURM_CPUS_ON_NODE ${tumor}
  """
}


process svcall {
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'
  input:
  set pair_id,file(ssbam),file(ssidx) from svbam
  output:
  file("${pair_id}.delly.vcf.gz") into dellyvcf
  file("${pair_id}.sssv.sv.vcf.gz") into svvcf
  file("${pair_id}.sv.vcf.gz") into svintvcf
  file("${pair_id}.sv.annot.txt") into svannot
  script:
  """
  source /etc/profile.d/modules.sh
  bash $baseDir/process_scripts/variants/svcalling.sh -b $ssbam -r $index_path -p $pair_id
  """
}

process gatkgvcf {
  errorStrategy 'ignore'
  //publishDir "$params.output", mode: 'copy'

  input:
  set pair_id,file(gbam),file(gidx) from gatkbam
  output:
  set pair_id,file("${pair_id}.gatk.g.vcf") into gvcf
  set pair_id,file("${pair_id}.gatk.vcf.gz") into gatkvcf

  script:
  """
  source /etc/profile.d/modules.sh
  module load gatk/3.5 python/2.7.x-anaconda bedtools/2.25.0 snpeff/4.2 vcftools/0.1.14 
  java -Djava.io.tmpdir=./ -Xmx32g -jar \$GATK_JAR -R ${reffa} -D ${dbsnp} -T HaplotypeCaller -stand_call_conf 30 -stand_emit_conf 10.0 -A FisherStrand -A QualByDepth -A VariantType -A DepthPerAlleleBySample -A HaplotypeScore -A AlleleBalance -variant_index_type LINEAR -variant_index_parameter 128000 --emitRefConfidence GVCF -I ${gbam} -o ${pair_id}.gatk.g.vcf -nct 2
  java -Djava.io.tmpdir=./ -Xmx32g -jar \$GATK_JAR -R ${reffa} -D ${dbsnp} -T GenotypeGVCFs -o gatk.vcf -nt 4 --variant ${pair_id}.gatk.g.vcf
  vcf-annotate -n --fill-type gatk.vcf | bcftools norm -c s -f ${reffa} -w 10 -O z -o ${pair_id}.gatk.vcf.gz -
  tabix ${pair_id}.gatk.vcf.gz
  """
}

process mpileup {
  errorStrategy 'ignore'
  //publishDir "$params.output", mode: 'copy'

  input:
  set pair_id,file(gbam),file(gidx) from sambam
  output:
  //file("${pair_id}.sampanel.vcf.gz") into samfilt
  set pair_id,file("${pair_id}.sam.vcf.gz") into samvcf
  script:
  """
  source /etc/profile.d/modules.sh
  module load python/2.7.x-anaconda samtools/intel/1.3 bedtools/2.25.0 bcftools/intel/1.3 snpeff/4.2 vcftools/0.1.14
  samtools mpileup -t 'AD,DP,INFO/AD' -ug -Q20 -C50 -f ${reffa} ${gbam} | bcftools call -vmO z -o ${pair_id}.sam.ori.vcf.gz
  vcf-sort ${pair_id}.sam.ori.vcf.gz | vcf-annotate -n --fill-type | bcftools norm -c s -f ${reffa} -w 10 -O z -o ${pair_id}.sam.vcf.gz -
  """
}
process hotspot {
  errorStrategy 'ignore'
  //publishDir "$params.output", mode: 'copy'

  input:
  set pair_id,file(gbam),file(gidx) from hsbam
  output:
  set pair_id,file("${pair_id}.hotspot.vcf.gz") into hsvcf
  when:
  params.cancer == "detect"
  script:
  """
  source /etc/profile.d/modules.sh
  module load python/2.7.x-anaconda samtools/intel/1.3 bedtools/2.25.0 bcftools/intel/1.3 snpeff/4.2 vcftools/0.1.14
  samtools mpileup -d 99999 -t 'AD,DP,INFO/AD' -uf ${reffa} ${gbam} > ${pair_id}.mpi
  bcftools filter -i "AD[1]/DP > 0.01" ${pair_id}.mpi | bcftools filter -i "DP > 50" | bcftools call -m -A |vcf-annotate -n --fill-type |  bcftools norm -c s -f /project/shared/bicf_workflow_ref/GRCh38/genome.fa -w 10 -O z -o ${pair_id}.lowfreq.vcf.gz -
  java -jar \$SNPEFF_HOME/SnpSift.jar annotate ${index_path}/cosmic.vcf.gz ${pair_id}.lowfreq.vcf.gz | java -jar \$SNPEFF_HOME/SnpSift.jar filter "(CNT[*] >0)" - |bgzip > ${pair_id}.hotspot.vcf.gz
  """
}
process speedseq {
  errorStrategy 'ignore'
  //publishDir "$params.output", mode: 'copy'

  input:
  set pair_id,file(gbam),file(gidx) from ssbam
  output:
  set pair_id,file("${pair_id}.ssvar.vcf.gz") into ssvcf

  script:
  """
  source /etc/profile.d/modules.sh
  module load python/2.7.x-anaconda samtools/intel/1.3 bedtools/2.25.0 bcftools/intel/1.3 snpeff/4.2 speedseq/20160506 vcftools/0.1.14
  speedseq var -t \$SLURM_CPUS_ON_NODE -o ssvar ${reffa} ${gbam}
  vcf-annotate -n --fill-type ssvar.vcf.gz| bcftools norm -c s -f ${reffa} -w 10 -O z -o ${pair_id}.ssvar.vcf.gz -
  """
}
process platypus {
  errorStrategy 'ignore'
  //publishDir "$params.output", mode: 'copy'

  input:
  set pair_id,file(gbam),file(gidx) from platbam

  output:
  //file("${pair_id}.platpanel.vcf.gz") into platfilt
  set pair_id,file("${pair_id}.platypus.vcf.gz") into platvcf

  script:
  """
  source /etc/profile.d/modules.sh
  module load python/2.7.x-anaconda bedtools/2.25.0 snpeff/4.2 platypus/gcc/0.8.1 bcftools/intel/1.3 samtools/intel/1.3 vcftools/0.1.14
  Platypus.py callVariants --minMapQual=10 --mergeClusteredVariants=1 --nCPU=\$SLURM_CPUS_ON_NODE --bamFiles=${gbam} --refFile=${reffa} --output=platypus.vcf
  vcf-sort platypus.vcf |vcf-annotate -n --fill-type -n |bgzip > platypus.vcf.gz
  tabix platypus.vcf.gz
  bcftools norm -c s -f ${reffa} -w 10 -O z -o ${pair_id}.platypus.vcf.gz platypus.vcf.gz
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
  //publishDir "$params.output", mode: 'copy'

  input:
  set fname,file(ss),file(gatk),file(sam),file(plat),file(hs) from vcflist
  
  output:
  set fname,file("${fname}.union.vcf.gz") into union
  script:
  if (params.cancer == "detect")
  """
  source /etc/profile.d/modules.sh
  module load gatk/3.5 python/2.7.x-anaconda bedtools/2.25.0 snpeff/4.2 bcftools/intel/1.3 samtools/intel/1.3 vcftools/0.1.14
  bedtools multiinter -i ${gatk} ${sam} ${ss} ${plat} ${hs} -names gatk sam ssvar platypus hotspot |cut -f 1,2,3,5 | bedtools sort -i stdin | bedtools merge -c 4 -o distinct >  ${fname}_integrate.bed
  bedtools intersect -header -v -a ${hs} -b ${sam} |bedtools intersect -header -v -a stdin -b ${gatk} | bedtools intersect -header -v -a stdin -b ${ss} |  bedtools intersect -header -v -a stdin -b ${plat} | bgzip > ${fname}.hotspot.nooverlap.vcf.gz
  vcf-sort ${sam} |bgzip > sam.vcf.gz 
  tabix ${fname}.hotspot.nooverlap.vcf.gz
  tabix ${gatk}
  tabix sam.vcf.gz
  tabix ${ss}
  tabix ${plat}
  java -Xmx32g -jar \$GATK_JAR -R ${reffa} -T CombineVariants --filteredrecordsmergetype KEEP_UNCONDITIONAL --variant:gatk ${gatk} --variant:sam sam.vcf.gz --variant:ssvar ${ss} --variant:platypus ${plat} --variant:hotspot ${fname}.hotspot.nooverlap.vcf.gz -genotypeMergeOptions PRIORITIZE -priority sam,ssvar,gatk,platypus,hotspot -o ${fname}.int.vcf
  perl $baseDir/scripts/uniform_integrated_vcf.pl ${fname}.int.vcf
  bgzip ${fname}_integrate.bed
  tabix ${fname}_integrate.bed.gz
  bgzip ${fname}.uniform.vcf
  tabix ${fname}.uniform.vcf.gz
  bcftools annotate -a ${fname}_integrate.bed.gz --columns CHROM,FROM,TO,CallSet -h ${index_path}/CallSet.header ${fname}.uniform.vcf.gz | bgzip > ${fname}.union.vcf.gz
  """
  else
  """
  source /etc/profile.d/modules.sh
  module load gatk/3.5 python/2.7.x-anaconda bedtools/2.25.0 snpeff/4.2 bcftools/intel/1.3 samtools/intel/1.3
  module load vcftools/0.1.14
  bedtools multiinter -i ${gatk} ${sam} ${ss} ${plat} -names gatk sam ssvar platypus |cut -f 1,2,3,5 | bedtools sort -i stdin | bedtools merge -c 4 -o distinct >  ${fname}_integrate.bed
  vcf-sort ${sam} |bgzip > sam.vcf.gz
  tabix ${fname}.hotspot.nooverlap.vcf.gz
  tabix ${gatk}
  tabix sam.vcf.gz
  tabix ${ss}
  tabix ${plat}
  java -Xmx32g -jar \$GATK_JAR -R ${reffa} -T CombineVariants --filteredrecordsmergetype KEEP_UNCONDITIONAL -genotypeMergeOptions PRIORITIZE --variant:gatk ${gatk} --variant:sam sam.vcf.gz --variant:ssvar ${ss} --variant:platypus ${plat} -priority sam,ssvar,gatk,platypus -o ${fname}.int.vcf
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

  script:
  """
  source /etc/profile.d/modules.sh
  module load python/2.7.x-anaconda bedtools/2.25.0 snpeff/4.2 bcftools/1.4.1
  tabix ${unionvcf}
  bcftools annotate -Oz -a ${index_path}/ExAC.vcf.gz -o ${fname}.exac.vcf.gz --columns CHROM,POS,AC_Het,AC_Hom,AC_Hemi,AC_Adj,AN_Adj,AC_POPMAX,AN_POPMAX,POPMAX ${unionvcf}
  tabix ${fname}.exac.vcf.gz 
  java -Xmx10g -jar \$SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c \$SNPEFF_HOME/snpEff.config ${snpeff_vers} ${fname}.exac.vcf.gz | java -jar \$SNPEFF_HOME/SnpSift.jar annotate -id ${index_path}/dbSnp.vcf.gz -  | java -jar \$SNPEFF_HOME/SnpSift.jar annotate -info CLNSIG,CLNDSDB,CLNDSDBID,CLNDBN,CLNREVSTAT,CLNACC ${index_path}/clinvar.vcf.gz - | java -jar \$SNPEFF_HOME/SnpSift.jar annotate -info CNT ${index_path}/cosmic.vcf.gz - | java -Xmx10g -jar \$SNPEFF_HOME/SnpSift.jar dbnsfp -v -db ${index_path}/dbNSFP.txt.gz - |bgzip > ${fname}.annot.vcf.gz
  tabix ${fname}.annot.vcf.gz
  """
}
