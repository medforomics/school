# ABRA2

This container has all of the code to run ABRA2

## Getting Started

These instructions will cover usage information and for the docker container 

### Prerequisities


In order to run this container you'll need docker installed.

* [Windows](https://docs.docker.com/windows/started)
* [OS X](https://docs.docker.com/mac/started/)
* [Linux](https://docs.docker.com/linux/started/)

### Usage

#### Container Parameters

To realign samples in a tumor/normal pair

abraopts are command line options for ABRA

```shell
opt=''
threads=`nproc`
ioopt="--in tumor.bam --out tumor.abra2.bam"
abraopts="--mbq 150 --mnf 5 --mer 0.05"

if [ -n "targetpanel.bed" ]
then		  
    opt="--targets targetpanel.bed"
fi
if [ -n "$nbam" ]
then
    ioopt="--in normal.bam,tumor.bam --out normal.abra2.bam,tumor.abra2.bam"
fi

for i in *.bam
do 
  docker run -v ${PWD}:/data docker.io/goalconsortium/abra2:latest samtools index -@ $threads $i
done

docker run -v ${PWD}:/data docker.io/goalconsortium/abra2:latest java -Xmx16G -jar /usr/local/bin/abra2.jar $ioopt --ref dnaref/genome.fa --threads $threads $opt --tmpdir tmpdir $abraopts > abra.log
```

You can run this container interactively using

```shell
docker run -it -v ${PWD}:/data docker.io/goalconsortium/abra2:latest bash
```

#### Environment Variables

- PATH "$PATH:/usr/local/bin"

#### Volumes

* `/data` - Working Directory

#### Useful File Locations

* `/usr/local/bin/abra2.jar` - ABRA JAR File
  
## Built With

* ABRA2 2.20
* SamTools 1.10

## Find Us

* [Docker](https://hub.docker.com/repository/docker/org/goalconsortium)

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the 
[tags on this repository](https://hub.docker.com/repository/docker/goalconsortium/abra/tags). 

## Authors

* **Jeremy Mathews
* **Brandi Cantarel, PhD

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

## Acknowledgments

* Funding for this work was provided by Cancer Prevention and Research Institute of Texas (RP150596).
