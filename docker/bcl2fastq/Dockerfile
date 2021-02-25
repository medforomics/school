FROM continuumio/miniconda3
MAINTAINER Brandi Cantarel

RUN conda install -c bih-cubi bcl2fastq2

#RUN wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/software/bcl2fastq/bcl2fastq2-v2-20-0-tar.zip; unzip bcl2fastq2-v2-20-0-tar.zip; tar xzvf bcl2fastq2-v2.20.0.422-Source.tar.gz; export C_INCLUDE_PATH=/usr/include/x86_64-linux-gnu; \ 
# ./bcl2fastq/src/configure --prefix=/usr/local/; \
# make; \
# make install; \
# rm -r /tmp/*; \

WORKDIR /
ENTRYPOINT ["bcl2fastq"]
CMD ["--version"]

