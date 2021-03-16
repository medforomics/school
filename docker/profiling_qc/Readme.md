# QC of VCF and BAM File and Variant Profiling

This container has all of the code to (i) determine quality metrics of BAM files including fastqc, mapping rate, on target percent and average coverage, (ii) compare tumor/normal pairs to ensure origins from the same patient, (iii) Microsattelite profiling and (iv) calculate the the bases per position.

## Getting Started

These instructions will cover usage information and for the docker container

### Prerequisities

In order to run this container you'll need docker installed.

* [Windows](https://docs.docker.com/windows/started)
* [OS X](https://docs.docker.com/mac/started/)
* [Linux](https://docs.docker.com/linux/started/)

### Usage

#### Container Parameters

sampleid is the Read Group/Sample Name
USER adds the user running the code to the QC stat file

To run DNA BAM QC
```shell
docker run -v ${PWD}:/data docker.io/goalconsortium/profiling_qc:1.0.0 bash /seqprg/school/process_scripts/alignment/bamqc.sh -c targetpanel.bed -n dna -r reference_files -b file.bam -p ${sampleid} -u $USER
```
To run RNA BAM QC
```shell
docker run -v ${PWD}:/data docker.io/goalconsortium/profiling_qc:1.0.0 bash /seqprg/school/process_scripts/alignment/bamqc.sh -p ${sampleid} -b file.bam -n rna 
```

To calculate the allele t each position in your BED file

```shell
docker run -v ${PWD}:/data docker.io/goalconsortium/profiling_qc:1.0.0 bam-readcount -l targetpanel.bed -w 0 -q 0 -b 25 -f genome.fa file.bam > file.bamreadcount.txt
```

Run NGS Checkmate to compare tumor/normal pair, determine if same patient.

dnaref is a directory with the genome sequence file (genome.fa)
caseid is the patientid, used for file naming
NGSCheckMate.bed -- can be downloaded from the NGSCheckMate program

```shell
docker run -v ${PWD}:/data docker.io/goalconsortium/vcfannot:0.5.36 bash /seqprg/school/process_scripts/variants/checkmate.sh -r dnaref -p ${caseid} -c dnaref/NGSCheckMate.bed -f
```
To run MSI-Sensor Pro

dnaref shoudl contain the microsattelite list file 

if Tumor/Normal pair add "-n normal.bam"

```shell
docker run -v ${PWD}:/data docker.io/goalconsortium/profiling_qc:1.0.0 bash /seqprg/school/process_scripts/variants/msisensor.sh -r dnaref -p ${caseid} -b tumor.bam -c targetpanel.bed
```

You can run this container interactively using

```shell
docker run -it -v ${PWD}:/data docker.io/goalconsortium/profiling_qc:1.0.0 bash
```

#### Environment Variables

- PATH "$PATH:/usr/local/bin/FastQC:/usr/local/bin"
- NCM_HOME "/usr/local/bin"
- PICARD "/usr/local/bin"
  - location of PICARD JAR File
- gitv "version_1.0.0"

#### Volumes

* `/data` - Working Directory

#### Useful File Locations

* `/usr/local/bin/picard.jar` - PICARD JAR File
* `/usr/local/bin/ncm.py` - ngs checkmate python script

## Built With

* Python 2.7.18
* SamTools + BCFTools + htslib 1.10
* bedtools 2.29.2
* fastqc_v0.11.5
* picard 2.21.7
* ngscheckmate v1.0.0
* msisensor-pro
* bamreadct v0.8.0

## Find Us

* [Docker](https://hub.docker.com/repository/docker/orgs/goalconsortium)
* Docker Container [goalconsortium/profiling_qc](https://hub.docker.com/repository/docker/goalconsortium/profiling_qc/general)
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
