sample=config["sample"]
fq1=config["fq1"]
fq2=config["fq2"]
indexPath="/refs"
captureBed="/panel/"+config["captureBed"]
itdBed="/panel/"+config["itdBed"]
if not config["mutectPon"] or config["mutectPon"].isspace():
  ponOpt=""
else:
  ponOpt="-q /panel/" + config["mutectPon"]
captureDir="/panel/"+config["captureDir"]
markdupsAlg=config["markdupsAlg"]
if markdupsAlg=="fgbio_umi":
  alignOpts="-u"
else:
  alignOpts=""
svAlgos=config["svAlgos"].split(' ')
snvIndelAlgos=config["snvIndelAlgos"].split(' ')
snpEff=config["snpEff"]
repoDir="/seqprg"
qcFiles=["_exoncoverage", ".genomecov", "_lowcoverage", ".covhist", ".flagstat", ".libcomplex", ".ontarget.flagstat"]

def getSVInput(wildcards):
  svInput = []
  if "cnvkit" in svAlgos:
    svInput += [
      "dnacallset/{sample}.call.cns".format(sample=sample),
      "dnacallset/{sample}.cns".format(sample=sample),
      "dnacallset/{sample}.cnr".format(sample=sample),
      "dnacallset/{sample}.answerplot.cnr".format(sample=sample),
      "dnacallset/{sample}.answerplot.cns".format(sample=sample),
      "dnacallset/{sample}.ballelefreq.txt".format(sample=sample),
      "dnacallset/{sample}.cnv.answer.txt".format(sample=sample),
      "dnacallset/{sample}.cnv.scatter.pdf".format(sample=sample),
    ]
  if "itdseek" in svAlgos:
    svInput += [
      "dnacallset/{sample}.itdseek_tandemdup.vcf.gz".format(sample=sample),
    ]
  if "svaba" in svAlgos:
    svInput += [
      "dnacallset/{sample}.svaba.vcf.gz".format(sample=sample),
      "dnacallset/{sample}.svaba.ori.vcf.gz".format(sample=sample),
      "dnacallset/{sample}.svaba.genefusion.txt".format(sample=sample),
    ]
  if "delly" in svAlgos:
    svInput += [
      "dnacallset/{sample}.delly.vcf.gz".format(sample=sample),
      "dnacallset/{sample}.delly.ori.vcf.gz".format(sample=sample),
      "dnacallset/{sample}.delly.genefusion.txt".format(sample=sample),
    ]
  return svInput

def getSNVInput(wildcards):
  snvInput = []
  if "pindel" in snvIndelAlgos:
    snvInput += [
      "dnacallset/{sample}.pindel_tandemdup.vcf.gz".format(sample=sample),
      "dnacallset/{sample}.pindel.sv.vcf.gz".format(sample=sample),
      "dnacallset/{sample}.pindel.genefusion.txt".format(sample=sample),
    ]
  # Other vcfs are handled by integrate rule
  return snvInput

rule all:
  input:
    "job.log",
    "{sample}.trim.R1.fastq.gz".format(sample=sample),
    "{sample}.trim.R2.fastq.gz".format(sample=sample),
    "{sample}.trimreport.txt".format(sample=sample),
    "dnaout/{sample}.bam".format(sample=sample),
    "dnaout/{sample}.bam.bai".format(sample=sample),
    "dnaout/{sample}.consensus.bam".format(sample=sample),
    "dnaout/{sample}.consensus.bam.bai".format(sample=sample),
    "{sample}.abra2.bam".format(sample=sample),
    "{sample}.abra2.bam.bai".format(sample=sample),
    expand("dnaout/{sample}{qcFile}.txt", sample=sample, qcFile=qcFiles),
    "dnacallset/{sample}.msi".format(sample=sample),
    "{sample}.final.bam".format(sample=sample),
    "{sample}.final.bam.bai".format(sample=sample),
    getSVInput,
    getSNVInput,
    "dnavcf/{sample}.union.vcf.gz".format(sample=sample),

rule init:
  output:
    "job.log"
  singularity:
    "docker://goalconsortium/trim_galore:1.1.3"
  shell:
    """
    mkdir -p dnaout
    mkdir -p dnacallset
    mkdir -p dnavcf
    touch job.log
    echo "Sample: {sample}" >> job.log
    echo "Read 1: {fq1}" >> job.log
    echo "Read 2: {fq2}" >> job.log
    echo "Capture Bed: {captureBed}" >> job.log
    echo "ITD Bed: {itdBed}" >> job.log
    echo "Mutect PoN Option: {ponOpt}" >> job.log
    echo "CNV Capture Dir: {captureDir}" >> job.log
    echo "MarkDups Alg: {markdupsAlg}" >> job.log
    echo "Align Opt: {alignOpts}" >> job.log
    echo "SV Algos: {svAlgos}" >> job.log
    echo "SNV/Indel Algos: {snvIndelAlgos}" >> job.log
    echo "snpEff Version: {snpEff}" >> job.log
    echo "SLURM_CPUS_ON_NODE: $SLURM_CPUS_ON_NODE" >> job.log
    """

rule dnaTrim:
  input:
    sequence1=expand("{fq1}", fq1=fq1),
    sequence2=expand("{fq2}", fq2=fq2),
    log="job.log",
  output:
    expand("{sample}.trim.R1.fastq.gz", sample=sample),
    expand("{sample}.trim.R2.fastq.gz", sample=sample),
    expand("{sample}.trimreport.txt", sample=sample),
    temp(expand("{fq1}_trimming_report.txt", fq1=fq1)),
    temp(expand("{fq2}_trimming_report.txt", fq2=fq2)),
  singularity:
    "docker://goalconsortium/trim_galore:1.1.3"
  shell:
    """
    echo $(date -u) "DNA Trimming Start" >> {input.log}
    bash {repoDir}/process_scripts/preproc_fastq/trimgalore.sh -f -p {sample} {input.sequence1} {input.sequence2}
    echo $(date -u) "DNA Trimming Finish" >> {input.log}
    """

rule dnaAlign:
  input:
    sequence1=expand("{sample}.trim.R1.fastq.gz", sample=sample),
    sequence2=expand("{sample}.trim.R2.fastq.gz", sample=sample),
  output:
    expand("dnaout/{sample}.bam", sample=sample),
    expand("dnaout/{sample}.bam.bai", sample=sample),
    temp(directory("dnaAlign_run")),
  singularity:
    "docker://goalconsortium/dna_alignment:1.0.9"
  shell:
    """
    echo $(date -u) "DNA Alignment Start" >> job.log
    mkdir dnaAlign_run
    cd dnaAlign_run
    ln -s ../{input.sequence1} {sample}.trim.R1.fastq.gz
    ln -s ../{input.sequence2} {sample}.trim.R2.fastq.gz
    bash {repoDir}/process_scripts/alignment/dnaseqalign.sh -r {indexPath} -p {sample} -x {input.sequence1} -y {input.sequence2} {alignOpts}
    mv {sample}.bam ../dnaout/{sample}.bam
    mv {sample}.bam.bai ../dnaout/{sample}.bam.bai

    echo $(date -u) "DNA Alignment End" >> ../job.log
    """

rule abra2:
  input:
    bam=expand("dnaout/{sample}.bam", sample=sample),
    bai=expand("dnaout/{sample}.bam.bai", sample=sample),
  output:
    expand("{sample}.abra2.bam", sample=sample),
    expand("{sample}.abra2.bam.bai", sample=sample),
    temp("abra.log"),
    temp(directory("tmpdir")),
  singularity:
    "docker://goalconsortium/abra2:1.0.9"
  shell:
    """
    echo $(date -u) "Abra2 Start" >> job.log
    bash {repoDir}/process_scripts/alignment/abra2.sh -r {indexPath} -p {sample} -b {input.bam} -c {captureBed}
    echo $(date -u) "Abra2 End" >> job.log
    """

rule markdups:
  input:
    bam=expand("dnaout/{sample}.bam", sample=sample),
    bai=expand("dnaout/{sample}.bam.bai", sample=sample),
  output:
    expand("dnaout/{sample}.consensus.bam", sample=sample),
    expand("dnaout/{sample}.consensus.bam.bai", sample=sample),
  singularity:
    "docker://goalconsortium/dna_alignment:1.0.9"
  shell:
    """
    echo $(date -u) "MarkDups Start" >> job.log
    bash {repoDir}/process_scripts/alignment/markdups.sh -a {markdupsAlg} -b {input.bam} -p {sample} -r {indexPath}
    mv {sample}.dedup.bam dnaout/{sample}.consensus.bam
    mv {sample}.dedup.bam.bai dnaout/{sample}.consensus.bam.bai
    echo $(date -u) "MarkDups End" >> job.log
    """

rule dna_bamqc:
  input:
    bam=expand("dnaout/{sample}.bam", sample=sample),
    bai=expand("dnaout/{sample}.bam.bai", sample=sample),
    trimreport=expand("{sample}.trimreport.txt", sample=sample)
  params:
    version="v4",
  output:
    #expand("dnaout/{sample}_fastqc.html", sample=sample),
    #expand("dnaout/{sample}_fastqc.zip", sample=sample),
    expand("dnaout/{sample}{qcFile}.txt", qcFile=qcFiles, sample=sample),
    temp(expand("{sample}.ontarget.bam", sample=sample)),
    temp(expand("{sample}.ontarget.bam.bai", sample=sample)),
    temp("panel.sorted.bed"),
  singularity:
    "docker://goalconsortium/profiling_qc:1.1.3"
  shell:
    """
    echo $(date -u) "QC Start" >> job.log
    bash {repoDir}/process_scripts/alignment/bamqc.sh -c {captureBed} -n dna -r {indexPath} -b {input.bam} -p {sample} -e {params.version}
    mv {sample}_exoncoverage.txt dnaout/{sample}_exoncoverage.txt
    mv {sample}.genomecov.txt dnaout/{sample}.genomecov.txt
    mv {sample}_lowcoverage.txt dnaout/{sample}_lowcoverage.txt
    mv {sample}.covhist.txt dnaout/{sample}.covhist.txt
    mv {sample}.flagstat.txt dnaout/{sample}.flagstat.txt
    mv {sample}.libcomplex.txt dnaout/{sample}.libcomplex.txt
    mv {sample}.ontarget.flagstat.txt dnaout/{sample}.ontarget.flagstat.txt
    #mv {sample}_fastqc.html dnaout/{sample}_fastqc.html
    #mv {sample}_fastqc.zip dnaout/{sample}_fastqc.zip

    echo $(date -u) "QC End" >> job.log
    """

rule msi:
  input:
    bam=expand("dnaout/{sample}.bam", sample=sample),
    bai=expand("dnaout/{sample}.bam.bai", sample=sample),
  output:
    expand("dnacallset/{sample}.msi", sample=sample),
    temp(expand("{sample}.msi_all", sample=sample)),
    temp(expand("{sample}.msi_dis", sample=sample)),
    temp(expand("{sample}.msi_unstable", sample=sample)),
  singularity:
    "docker://goalconsortium/profiling_qc:1.1.3"
  shell:
    """
    echo $(date -u) "MSI Start" >> job.log
    bash {repoDir}/process_scripts/variants/msisensor.sh -r {indexPath} -p {sample} -b {input.bam} -c {captureBed}
    mv {sample}.msi dnacallset/{sample}.msi
    echo $(date -u) "MSI End" >> job.log
    """

if "itdseek" in svAlgos:
  rule itdseek:
    input:
      bam=expand("dnaout/{sample}.bam", sample=sample),
      bai=expand("dnaout/{sample}.bam.bai", sample=sample),
    output:
      expand("dnacallset/{sample}.itdseek_tandemdup.vcf.gz", sample=sample),
      temp(directory("itdseek_run")),
    singularity:
      "docker://goalconsortium/structuralvariant:1.1.4"
    shell:
      """
      echo $(date -u) "ITDSeek Start" >> job.log
      mkdir itdseek_run
      cd itdseek_run
      ln -s ../{input.bam} {sample}.bam
      ln -s ../{input.bai} {sample}.bam.bai
      bash {repoDir}/process_scripts/variants/svcalling.sh -b {sample}.bam -r {indexPath} -p {sample} -l {itdBed} -a itdseek -g {snpEff} -f
      mv {sample}.itdseek_tandemdup.vcf.gz ../dnacallset/{sample}.itdseek_tandemdup.vcf.gz

      echo $(date -u) "ITDSeek End" >> ../job.log
      """

if "delly" in svAlgos:
  rule delly:
    input:
      bam=expand("dnaout/{sample}.bam", sample=sample),
      bai=expand("dnaout/{sample}.bam.bai", sample=sample),
    output:
      expand("dnacallset/{sample}.delly.vcf.gz", sample=sample),
      expand("dnacallset/{sample}.delly.ori.vcf.gz", sample=sample),
      expand("dnacallset/{sample}.delly.genefusion.txt", sample=sample),
      temp(directory("delly_run")),
    singularity:
      "docker://goalconsortium/structuralvariant:1.1.4"
    shell:
      """
        echo $(date -u) "Delly Start" >> job.log
        mkdir delly_run
        cd delly_run
        ln -s ../{input.bam} {sample}.bam
        ln -s ../{input.bai} {sample}.bam.bai
        bash {repoDir}/process_scripts/variants/svcalling.sh -r {indexPath} -b {sample}.bam -p {sample} -a delly -g {snpEff} -f
        mv {sample}.delly.vcf.gz ../dnacallset/{sample}.delly.vcf.gz
        mv {sample}.delly.ori.vcf.gz ../dnacallset/{sample}.delly.ori.vcf.gz
        mv {sample}.delly.genefusion.txt ../dnacallset/{sample}.delly.genefusion.txt
        echo $(date -u) "Delly End" >> ../job.log
      """

if "svaba" in svAlgos:
  rule svaba:
    input:
      bam=expand("dnaout/{sample}.bam", sample=sample),
      bai=expand("dnaout/{sample}.bam.bai", sample=sample),
      consensusBam=expand("dnaout/{sample}.consensus.bam", sample=sample),
      consensusBai=expand("dnaout/{sample}.consensus.bam.bai", sample=sample)
    output:
      expand("dnacallset/{sample}.svaba.vcf.gz", sample=sample),
      expand("dnacallset/{sample}.svaba.ori.vcf.gz", sample=sample),
      expand("dnacallset/{sample}.svaba.genefusion.txt", sample=sample),
      temp(directory("svaba_run")),
    singularity:
      "docker://goalconsortium/structuralvariant:1.1.4"
    shell:
      """
        echo $(date -u) "Svaba Start" >> job.log
        mkdir svaba_run
        cd svaba_run
        ln -s ../{input.bam} {sample}.bam
        ln -s ../{input.bai} {sample}.bam.bai
        ln -s ../{input.consensusBam} {sample}.consensus.bam
        ln -s ../{input.consensusBai} {sample}.consensus.bam.bai
        bash {repoDir}/process_scripts/variants/svcalling.sh -r {indexPath} -b {sample}.bam -p {sample} -a svaba -g {snpEff} -f
        mv {sample}.svaba.vcf.gz ../dnacallset/{sample}.svaba.vcf.gz
        mv {sample}.svaba.ori.vcf.gz ../dnacallset/{sample}.svaba.ori.vcf.gz
        mv {sample}.svaba.genefusion.txt ../dnacallset/{sample}.svaba.genefusion.txt

        echo $(date -u) "Svaba End" >> ../job.log
      """

if "cnvkit" in svAlgos:
  rule cnvkit:
    input:
      bam=expand("dnaout/{sample}.bam", sample=sample),
      bai=expand("dnaout/{sample}.bam.bai", sample=sample)
    output:
      expand("dnacallset/{sample}.call.cns", sample=sample),
      expand("dnacallset/{sample}.cns", sample=sample),
      expand("dnacallset/{sample}.cnr", sample=sample),
      expand("dnacallset/{sample}.answerplot.cnr", sample=sample),
      expand("dnacallset/{sample}.answerplot.cns", sample=sample),
      expand("dnacallset/{sample}.ballelefreq.txt", sample=sample),
      expand("dnacallset/{sample}.cnv.answer.txt", sample=sample),
      expand("dnacallset/{sample}.cnv.scatter.pdf", sample=sample),
      temp(directory("cnv_run")),
    singularity:
      "docker://goalconsortium/structuralvariant:1.1.4"
    shell:
      """
        echo $(date -u) "CNV Start" >> job.log
        mkdir cnv_run
        cd cnv_run
        ln -s ../{input.bam} {sample}.bam
        ln -s ../{input.bai} {sample}.bam.bai
        bash {repoDir}/process_scripts/variants/cnvkit.sh -r {indexPath} -b {sample}.bam  -p {sample} -d {captureDir}
        mv {sample}.call.cns ../dnacallset/{sample}.call.cns
        mv {sample}.cns ../dnacallset/{sample}.cns
        mv {sample}.cnr ../dnacallset/{sample}.cnr
        mv {sample}.answerplot.cnr ../dnacallset/{sample}.answerplot.cnr
        mv {sample}.answerplot.cns ../dnacallset/{sample}.answerplot.cns
        mv {sample}.ballelefreq.txt ../dnacallset/{sample}.ballelefreq.txt
        mv {sample}.cnv.answer.txt ../dnacallset/{sample}.cnv.answer.txt
        mv {sample}.cnv.scatter.pdf ../dnacallset/{sample}.cnv.scatter.pdf
        echo $(date -u) "CNV End" >> ../job.log
      """

if "pindel" in snvIndelAlgos:
  rule pindel:
    input:
      bam=expand("dnaout/{sample}.consensus.bam", sample=sample),
      bai=expand("dnaout/{sample}.consensus.bam.bai", sample=sample),
    output:
      expand("dnacallset/{sample}.pindel_tandemdup.vcf.gz", sample=sample),
      expand("dnacallset/{sample}.pindel.vcf.gz", sample=sample),
      expand("dnacallset/{sample}.pindel.sv.vcf.gz", sample=sample),
      expand("dnacallset/{sample}.pindel.genefusion.txt", sample=sample),
      temp(directory("pindel_run")),
    singularity:
      "docker://goalconsortium/structuralvariant:1.1.4"
    shell:
      """
      echo $(date -u) "Pindel Start" >> job.log
      mkdir pindel_run
      cd pindel_run
      ln -s ../{input.bam} {sample}.consensus.bam
      ln -s ../{input.bai} {sample}.consensus.bam.bai
      bash {repoDir}/process_scripts/variants/svcalling.sh -r {indexPath} -p {sample} -l {itdBed} -a pindel -c {captureBed} -g {snpEff} -f
      mv {sample}.pindel_tandemdup.vcf.gz ../dnacallset/{sample}.pindel_tandemdup.vcf.gz
      mv {sample}.pindel.vcf.gz ../dnacallset/{sample}.pindel.vcf.gz
      mv {sample}.pindel.sv.vcf.gz ../dnacallset/{sample}.pindel.sv.vcf.gz
      mv {sample}.pindel.genefusion.txt ../dnacallset/{sample}.pindel.genefusion.txt
      echo $(date -u) "Pindel End" >> ../job.log
      """

if "mutect" in snvIndelAlgos:
  rule gatkbam:
    input:
      bam=expand("dnaout/{sample}.consensus.bam", sample=sample),
      bai=expand("dnaout/{sample}.consensus.bam.bai", sample=sample),
    output:
      expand("{sample}.final.bam", sample=sample),
      expand("{sample}.final.bam.bai", sample=sample),
      temp(directory("gatkbam_run")),
    singularity:
      "docker://goalconsortium/variantcalling:1.1.5"
    shell:
      """
      echo $(date -u) "GATKBam Start" >> job.log
      mkdir gatkbam_run
      cd gatkbam_run
      ln -s ../{input.bam} {sample}.consensus.bam
      ln -s ../{input.bai} {sample}.consensus.bam.bai
      bash {repoDir}/process_scripts/variants/gatkrunner.sh -a gatkbam -b {sample}.consensus.bam -r {indexPath} -p {sample}
      mv {sample}.final.bam ../{sample}.final.bam
      mv {sample}.final.bam.bai ../{sample}.final.bam.bai
      echo $(date -u) "GATKBam End" >> ../job.log
      """

  rule mutect:
    input:
      bam=expand("{sample}.final.bam", sample=sample),
      bai=expand("{sample}.final.bam.bai", sample=sample),
    output:
      expand("dnacallset/{sample}.mutect.vcf.gz", sample=sample),
      expand("dnacallset/{sample}.mutect.ori.vcf.gz", sample=sample),
      temp(directory("mutect_run")),
    singularity:
      "docker://goalconsortium/variantcalling:1.1.5"
    shell:
      """
        echo $(date -u) "Mutect Start" >> job.log
        mkdir mutect_run
        cd mutect_run
        ln -s ../{input.bam} {sample}.final.bam
        ln -s ../{input.bai} {sample}.final.bam.bai
        bash {repoDir}/process_scripts/variants/germline_vc.sh {ponOpt} -r {indexPath} -p {sample} -b {captureBed} -a mutect
        bash {repoDir}/process_scripts/variants/uni_norm_annot.sh -g {snpEff} -r {indexPath} -p {sample}.mutect -v {sample}.mutect.vcf.gz

        mv {sample}.mutect.vcf.gz ../dnacallset/{sample}.mutect.vcf.gz
        mv {sample}.mutect.ori.vcf.gz ../dnacallset/{sample}.mutect.ori.vcf.gz

        echo $(date -u) "Mutect End" >> ../job.log
      """

if "strelka2" in snvIndelAlgos:
  rule germstrelka:
    input:
      bam=expand("dnaout/{sample}.consensus.bam", sample=sample),
      bai=expand("dnaout/{sample}.consensus.bam.bai", sample=sample),
    output:
      expand("dnacallset/{sample}.strelka2.vcf.gz", sample=sample),
      expand("dnacallset/{sample}.strelka2.ori.vcf.gz", sample=sample),
      temp(directory("strelka_run")),
    singularity:
      "docker://goalconsortium/variantcalling:1.1.5"
    shell:
      """
        echo $(date -u) "Strelka Start" >> job.log
        mkdir strelka_run
        cd strelka_run
        ln -s ../{input.bam} {sample}.consensus.bam
        ln -s ../{input.bai} {sample}.consensus.bam.bai
        bash {repoDir}/process_scripts/variants/germline_vc.sh -r {indexPath} -p {sample} -a strelka2 -b {captureBed}
        bash {repoDir}/process_scripts/variants/uni_norm_annot.sh -g {snpEff} -r {indexPath} -p {sample}.strelka2 -v {sample}.strelka2.vcf.gz

        mv {sample}.strelka2.vcf.gz ../dnacallset/{sample}.strelka2.vcf.gz
        mv {sample}.strelka2.ori.vcf.gz ../dnacallset/{sample}.strelka2.ori.vcf.gz

        echo $(date -u) "Strelka End" >> ../job.log
      """

if "platypus" in snvIndelAlgos:
  rule germplatypus:
    input:
      bam=expand("dnaout/{sample}.consensus.bam", sample=sample),
      bai=expand("dnaout/{sample}.consensus.bam.bai", sample=sample),
    output:
      expand("dnacallset/{sample}.platypus.vcf.gz", sample=sample),
      expand("dnacallset/{sample}.platypus.ori.vcf.gz", sample=sample),
      temp(directory("platypus_run")),
    singularity:
      "docker://goalconsortium/variantcalling:1.1.5"
    shell:
      """
        echo $(date -u) "Platypus Start" >> job.log
        mkdir platypus_run
        cd platypus_run
        ln -s ../{input.bam} {sample}.consensus.bam
        ln -s ../{input.bai} {sample}.consensus.bam.bai
        bash {repoDir}/process_scripts/variants/germline_vc.sh -r {indexPath} -p {sample} -a platypus -b {captureBed}
        bash {repoDir}/process_scripts/variants/uni_norm_annot.sh -g {snpEff} -r {indexPath} -p {sample}.platypus -v {sample}.platypus.vcf.gz

        mv {sample}.platypus.vcf.gz ../dnacallset/{sample}.platypus.vcf.gz
        mv {sample}.platypus.ori.vcf.gz ../dnacallset/{sample}.platypus.ori.vcf.gz
        echo $(date -u) "Platypus End" >> ../job.log
      """

if "fb" in snvIndelAlgos:
  rule germfb:
    input:
      bam=expand("dnaout/{sample}.consensus.bam", sample=sample),
      bai=expand("dnaout/{sample}.consensus.bam.bai", sample=sample),
    output:
      expand("dnacallset/{sample}.fb.vcf.gz", sample=sample),
      expand("dnacallset/{sample}.fb.ori.vcf.gz", sample=sample),
      temp(directory("fb_run")),
    singularity:
      "docker://goalconsortium/variantcalling:1.1.5"
    shell:
      """
        echo $(date -u) "FreeBayes Start" >> job.log
        mkdir fb_run
        cd fb_run
        ln -s ../{input.bam} {sample}.consensus.bam
        ln -s ../{input.bai} {sample}.consensus.bam.bai
        bash {repoDir}/process_scripts/variants/germline_vc.sh -r {indexPath} -p {sample} -a fb -b {captureBed}
        bash {repoDir}/process_scripts/variants/uni_norm_annot.sh -g {snpEff} -r {indexPath} -p {sample}.fb -v {sample}.fb.vcf.gz

        mv {sample}.fb.vcf.gz ../dnacallset/{sample}.fb.vcf.gz
        mv {sample}.fb.ori.vcf.gz ../dnacallset/{sample}.fb.ori.vcf.gz
        echo $(date -u) "FreeBayes End" >> ../job.log
      """

rule integrate:
  input:
    vcfs=expand("dnacallset/{sample}.{snvIndelAlgos}.vcf.gz", sample=sample, snvIndelAlgos=snvIndelAlgos),
  output:
    expand("dnavcf/{sample}.union.vcf.gz", sample=sample),
    temp(directory("integrate_run")),
  singularity:
    "docker://goalconsortium/structuralvariant:1.1.4"
  shell:
    """
      echo $(date -u) "Integrate Start" >> job.log
      mkdir integrate_run
      cd integrate_run
      for algo in {input.vcfs}
      do
          ln -s ../$algo $(basename $algo)
      done
      bash {repoDir}/process_scripts/variants/union.sh -r {indexPath} -p {sample}
      mv {sample}.union.vcf.gz ../dnavcf/{sample}.union.vcf.gz
      echo $(date -u) "Integrate End" >> ../job.log
    """
