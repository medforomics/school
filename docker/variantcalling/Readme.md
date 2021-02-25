# QC of VCF and BAM File and Variant Profiling

This container has all of the code to detect and annotate small SNVs (single nucleotide variants) and InDels (< 50 bp insertions/deletions)

## Getting Started

These instructions will cover usage information and for the docker container

### Prerequisities

In order to run this container you'll need docker installed.

* [Windows](https://docs.docker.com/windows/started)
* [OS X](https://docs.docker.com/mac/started/)
* [Linux](https://docs.docker.com/linux/started/)

### Usage

#### Container Parameters

dnaref -- contains dbSNP and genome.fa files
sampleid -- output file prefix

* To run GATK BQSR
```shell
docker run -v ${PWD}:/data docker.io/goalconsortium/variantcalling:1.0.0 bash /seqprg/school/process_scripts/variants/gatkrunner.sh -a gatkbam -b file.bam -r dnaref -p ${sampleid}
```
* To run small variant calling

Variant Calling Methods:
- fb: FreeBayes
- platypus: Platypus
- strelka2: Strelka2
- mutect: MuTect2 (GATK4)
- shimmer: Shimmer

dnaref -- contains dbSNP and genome.fa files

In Single Sample Mode

```shell
docker run -v ${PWD}:/data docker.io/goalconsortium/variantcalling:1.0.0 bash /seqprg/school/process_scripts/variants/germline_vc.sh -r dnaref -p $prefix -a $method -b targetpanel.bed
```

In Tumor/Normal Mode
A panel of normal VCF can be used by MUTECT to remove common artifacts 

```shell
docker run -v ${PWD}:/data docker.io/goalconsortium/variantcalling:1.0.0 bash /seqprg/school/process_scripts/variants/somatic_vc.sh -p $prefix -a $method -b targetpanel.bed -t tumor.bam -n normal.bam -q mutect.pon.vcf.gz
```

* To Annotate VCF Files

snpeff_vers - SNP Eff Genome Version Name
dnaref has genome.fa, genome.dict, annotation files including:
  - gnomad.txt.gz
  - oncokb_hotspot.txt.gz
  - repeat_regions.bed.gz
  - dbSnp.vcf.gz
  - clinvar.vcf.gz
  - cosmic.vcf.gz
  - dbNSFP.txt.gz

```shell
docker run -v ${PWD}:/data docker.io/goalconsortium/variantcalling:1.0.0 bash /seqprg/school/process_scripts/variants/uni_norm_annot.sh -g $snpeff_vers -r dnaref -p $prefix -v file.vcf.gz
```

* Union VCF -- needs a genome dictorary file called genome.dict

```shell
docker run -v ${PWD}:/data docker.io/goalconsortium/variantcalling:1.0.0 bash /seqprg/school/process_scripts/variants/union.sh -r reference directory -p ${caseid}
```

You can run this container interactively using

```shell
docker run -it -v ${PWD}:/data docker.io/goalconsortium/variantcalling:1.0.0 bash
```

#### Environment Variables

- PATH "$PATH:/usr/local/bin/strelka-2.9.10.centos6_x86_64/bin:/usr/local/bin"
- SNPEFF_HOME "/usr/local/bin/snpEff"
- PICARD "/usr/local/bin"
  - location of PICARD JAR File

#### Volumes

* `/data` - Working Directory

#### Useful File Locations

* `/usr/local/bin/picard.jar` - PICARD JAR File


## Built With

* python 2.7.18
* SamTools + BCFTools + htslib 1.10
* Platypus 0.8.1
* bedtools 2.29.2
* vcftools 0.1.14
* snpEff v4_3t
* picard 2.21.7
* strelka 2.9.10
* manta 1.6.0
* freebayes v1.2.0
* shimmer v0.1.1
* VT 0.57721

## Find Us

* [Docker](https://hub.docker.com/repository/docker/orgs/goalconsortium)
* Docker Container [goalconsortium/variantcalling](https://hub.docker.com/repository/docker/goalconsortium/variantcalling/general)
* Git Repo [SCHOOL](https://github.com/bcantarel/school)

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the 
[tags on this repository](https://hub.docker.com/repository/docker/goalconsortium/abra/tags). 

## Authors

* Jeremy Mathews
* Brandi Cantarel, PhD

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

## Acknowledgments

* Funding for this work was provided by Cancer Prevention and Research Institute of Texas (RP150596).
