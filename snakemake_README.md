## Introduction
This pipeline implements the tumor only version of the GOAL Consortium Consensus Pipeline in Snakemake. Documentation for the pipeline can be found at the following locations:

* [**Docker Containers**](https://hub.docker.com/orgs/goalconsortium)
* [**Pipeline Command Scripts (inside containers)**](https://github.com/medforomics/process_scripts)
* [**Nextflow and Postprocessing Scripts and Docker files**](https://github.com/GoalConsortium/school)
* [**DNA Nexus Applets**](https://github.com/GoalConsortium/dnanexus_applets)


## SNV/Indel Calling
The pipeline can be configured to report SNV and small indels calls from the following callers:

* [**Mutect2**](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2)
* [**Pindel**](http://gmt.genome.wustl.edu/packages/pindel/index.html)
* [**Platypus**](https://www.well.ox.ac.uk/research/research-groups/lunter-group/lunter-group/platypus-a-haplotype-based-variant-caller-for-next-generation-sequence-data)
* [**Strelka2**](https://github.com/Illumina/strelka)
* [**freebayes**](https://github.com/freebayes/freebayes)

## SV/CNV Calling
The pipeline can be configured to report structural variant and copy number calls from the following callers:

* [**Pindel**](http://gmt.genome.wustl.edu/packages/pindel/index.html)
* [**ITDseek**](https://github.com/tommyau/itdseek)
* [**Delly**](https://github.com/dellytools/delly)
* [**SvABA**](https://github.com/walaj/svaba)
* [**CNVkit**](https://cnvkit.readthedocs.io/en/stable/)

## Example Run

    export SINGULARITYENV_SLURM_CPUS_ON_NODE=$SLURM_CPUS_ON_NODE
    snakemake -r --cores all --use-singularity \
        --singularity-prefix /data/user/$USER/GOAL/.snakemake \
        --singularity-args "-B /path/to/panel/folder:/panel,/path/to/references/folder:/refs --cleanenv -c " \
        --directory /host/directory/to/work/in \
        --configfile /path/to/configFile \
        --config sample={Sample name} fq1={read one file} fq2={read two file}
        --config sample=$sample fq1=$fq1_file fq2=$fq2_file

## Configuring the Pipeline
### Environment
Singularity loads environment variables from the host to the container by default. This can create unintended side effects. To circumvent this, you should add --cleanenv flag to snakemake's singularity-args

    snakemake -r --cores ... --use-singularity \
        ... \
        --singularity-args "... --cleanenv ...

### Job Scheduler
An environment variable in the container to indicate the number of cores available to the job is necessary to parallelize computation. This can be set by running the following command in your job script prior to calling the pipelime. (If you are not using Slurm, then $SLURM_CPUS_ON_NODE must be substituted with your job scheduler's analog or a hard coded value)

    export SINGULARITYENV_SLURM_CPUS_ON_NODE=$SLURM_CPUS_ON_NODE

### SnpEff / Using GRCh37
The pipeline uses SnpEff to annotate variants. By default, the containers have the database for version GRCh38.86. If you want to use a different version, then you must download the desired database from http://pcingola.github.io/SnpEff/download/ and mount the database folder into the singularity containers.

    snakemake -r --cores ... --use-singularity \
        ... \
        --singularity-args "... -B /path/to/snpEff/data/GRCh37.75:/usr/local/bin/snpEff/data/GRCh37.75 ..."

### Reference Folder
The guide to creating the necessary reference files for each container can be found under each image repository on https://hub.docker.com/u/goalconsortium

* https://hub.docker.com/r/goalconsortium/trim_galore
* https://hub.docker.com/r/goalconsortium/dna_alignment
* https://hub.docker.com/r/goalconsortium/profiling_qc
* https://hub.docker.com/r/goalconsortium/structuralvariant
* https://hub.docker.com/r/goalconsortium/variantcalling
* https://hub.docker.com/r/goalconsortium/abra2

The references folder must be mounted to `/refs` inside the containers

    snakemake -r --cores all --use-singularity \
        ... \
        --singularity-args "... -B /path/to/references/folder:/refs ..."

### Panel Folder
A folder with all the panel bed files must be mounted into the containers. Note the filenames inside the captureDir are fixed. The other filenames can be specified in the pipeline config file.

    ./
        {captureBed}
        {itdBed}
        {mutectPon}
        {captureDir}/
            cnvkit.antitargets.bed
            cnvkit.targets.bed
            commonsnps.bed
            pon.cnn

The panel folder must be mounted to `/panel` inside the containers

    snakemake -r --cores all --use-singularity \
        ... \
        --singularity-args "... -B /path/to/panel/folder:/panel ..."


### Config file
The pipeline requires a config file to be set with the following keys.

* **captureBed:** The bed file in the panel folder for SNV/Indel calling with Mutect2, Pindel, Strelka2, Platypus, and freebayes
* **itdBed:** The bed file for ITD calling with Pindel and ITDseek.
* **mutectPon:** The Panel of Normal vcf to be used with Mutect2.
* **captureDir:** The name of the directory inside the panel folder with the Panel of Normal files to be used with CNVkit.
* **markdupsAlg:** The tool to be used to mark duplicates on a BAM. Options are **picard**, **picard_umi**, **samtools**, **fgbio_umi**.
* **svAlgos:** Space separated list of tools to be used to call structural variants. Options are **delly**, **svaba**, **cnvkit**, and **itdseek**.
* **snvIndelAlgos:** Space separated list of tools to be used to call SNVs and Indels variants. Options are **mutect**, **pindel**, **strelka2**, **platypus**, and **fb**
* **snpEff:** SnpEff version to be used to annotate.


*/configFile.json*

    {
      "captureBed": "goal_core497.hg19.bed",
      "itdBed": "itd_genes.bed",
      "mutectPon": "mutect2.pon.vcf.gz",
      "captureDir": "capture_dir",
      "markdupsAlg": "picard",
      "svAlgos": "delly svaba cnvkit itdseek",
      "snvIndelAlgos": "mutect pindel strelka2 platypus fb",
      "snpEff": "GRCh38.86"
    }

On CLI:

    snakemake -r --cores all \
        ... \
        --configfile /path/to/configFile

### Convention
When running the pipeline, the paired end FASTQ files should be copied, linked, or moved to a working directory for the analysis. The directory path and filenames for the read files must be passed as parameters to the pipeline. A sample name must be passed to the pipeline too for naming output files.

    snakemake -r --cores all \
        ... \
        --directory /host/directory/to/work/in \
        --config sample={Sample name} fq1={read one file in the directory} fq2={read one file in the directory}
