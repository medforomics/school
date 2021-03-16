# RNA_Alignment

This container has all of the code to create DNA alignments.  Can also sort, index and mark duplicates in alignments. 

## Getting Started

These instructions will cover usage information and for the docker container

### Prerequisities


In order to run this container you'll need docker installed.

* [Windows](https://docs.docker.com/windows/started)
* [OS X](https://docs.docker.com/mac/started/)
* [Linux](https://docs.docker.com/linux/started/)

### Usage

#### Container Parameters

To create a sorted BAM with HiSAT2 and samtools

rnaref is a directory with HiSAT2 indexed genome file
sampleid is the read group for the sample


```shell
docker run -v ${PWD}:/data docker.io/goalconsortium/dna_alignment:1.0.0 bash /seqprg/school/process_scripts/alignment/dnaseqalign.sh -r $humanref -p ${sampleid} -x seq.R1.fastq.gz -y seq.R2.fastq.gz
```

If you have a UMI sequence in the read name, you can transfer the UMI to the RX tag in the BAM using the -u option

```shell
docker run -v ${PWD}:/data docker.io/goalconsortium/dna_alignment:1.0.0 bash /seqprg/school/process_scripts/alignment/dnaseqalign.sh -r $humanref -p ${sampleid} -x seq.R1.fastq.gz -y seq.R2.fastq.gz -u 
```

You can run this container interactively using

```shell
docker run -it -v ${PWD}:/data docker.io/goalconsortium/rna_alignment:latest bash
```

#### Volumes

* `/data` - Working Directory

## Built With

* HiSAT 2.1.0
* SamTools 1.10
* python 2.7.18

## Find Us

* [Docker](https://hub.docker.com/repository/docker/orgs/goalconsortium)
* Docker Container [goalconsortium/rna_alignment](https://hub.docker.com/repository/docker/goalconsortium/rna_alignment/general)
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
