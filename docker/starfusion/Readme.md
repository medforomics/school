# Star-Fusion

This container can detect gene fusions in RNASeq data using Star-Fusion.

## Getting Started

These instructions will cover usage information and for the docker container

### Prerequisities

In order to run this container you'll need docker installed.

* [Windows](https://docs.docker.com/windows/started)
* [OS X](https://docs.docker.com/mac/started/)
* [Linux](https://docs.docker.com/linux/started/)

### Usage

#### Container Parameters

To Run Star Fusion and Annotate Gene Fusion Events

sampleid -- Sample Name
rnaref - [CTAT Resource Library Directory](https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/) (untarred unzip)

```shell
docker run -v ${PWD}:/data docker.io/goalconsortium/starfusion:1.0.0 bash /seqprg/school/process_scripts/alignment/starfusion.sh -p ${sampleid} -r $rnaref -a seq.R1.fastq.gz -b seq.trim.R2.fastq.gz -f 
```

You can run this container interactively using

```shell
docker run -it -v ${PWD}:/data docker.io/goalconsortium/starfusion:latest bash
```

#### Environment Variables

- PATH "$PATH:/usr/local/bin"
  
#### Volumes

* `/data` - Working Directory
  
## Built With

* Starfusion 1.9.0
* bedtools 2.29.2
* AGFusion

## Find Us

* [Docker](https://hub.docker.com/repository/docker/orgs/goalconsortium)
* Docker Container [goalconsortium/alignment](https://hub.docker.com/repository/docker/goalconsortium/starfusion/general)
* Git Repo [SCHOOL](https://github.com/bcantarel/school)

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the 
[tags on this repository](https://hub.docker.com/repository/docker/goalconsortium/starfusion/tags). 

## Authors

* Jeremy Mathews
* Brandi Cantarel, PhD

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

## Acknowledgments

* Funding for this work was provided by Cancer Prevention and Research Institute of Texas (RP150596).
