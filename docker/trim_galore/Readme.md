# Trim Galore

This container can trim reads.

## Getting Started

These instructions will cover usage information and for the docker container

### Prerequisities

In order to run this container you'll need docker installed.

* [Windows](https://docs.docker.com/windows/started)
* [OS X](https://docs.docker.com/mac/started/)
* [Linux](https://docs.docker.com/linux/started/)

### Usage

#### Container Parameters

To trim reads

```shell
docker run -v ${PWD}:/data docker.io/goalconsortium/trim_galore:1.0.0 -f -p ${sampleid} seq.R1.fastq.gz seq.R2.fastq.gz
```

You can run this container interactively using

```shell
docker run -it -v ${PWD}:/data docker.io/goalconsortium/trim_galore:latest bash
```

#### Environment Variables

- PATH "$PATH:/usr/local/bin"
  
#### Volumes

* `/data` - Working Directory

#### Useful File Locations

* `/usr/local/bin/picard.jar` - PICARD JAR File
  
## Built With

* trim_galore v0.4.1
* cutadapt 1.9.1

## Find Us

* [Docker](https://hub.docker.com/repository/docker/orgs/goalconsortium)
* Docker Container [goalconsortium/alignment](https://hub.docker.com/repository/docker/goalconsortium/trim_galore/general)
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
