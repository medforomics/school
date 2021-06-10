FROM ubuntu:16.04
LABEL Author="Jeremy Mathews" \
      Maintainer="genome.school@gmail.com"

# install additional requirements
RUN apt-get update; apt-get install -y build-essential checkinstall; \
    apt-get install -y libreadline-gplv2-dev libncursesw5-dev libssl-dev libsqlite3-dev tk-dev libgdbm-dev libc6-dev libbz2-dev zlib1g-dev liblzma-dev libcurl4-gnutls-dev libncurses5-dev wget unzip git default-jre default-jdk r-base cmake patch gfortran hdf5-tools libboost-date-time-dev libboost-program-options-dev libboost-system-dev libboost-filesystem-dev libboost-iostreams-dev libhdf5-dev fort77 xorg-dev libblas-dev gcc-multilib gobjc++ aptitude texlive-latex-base libcairo2-dev pigz;

# Make destination and work directory
RUN mkdir -p /opt;
RUN mkdir -p /opt/bin;
#ENV http_proxy http://proxy.swmed.edu:3128/
#ENV https_proxy http://proxy.swmed.edu:3128/

# Install Python
RUN cd /usr/src; wget https://www.python.org/ftp/python/2.7.18/Python-2.7.18.tgz; \
    tar -zxvf Python-2.7.18.tgz; rm Python-2.7.18.tgz; \
    cd Python-2.7.18; \
    ./configure --enable-optimizations; \
    make altinstall;
RUN cd /usr/local/bin; \
    mv python2.7 python;
RUN python -m ensurepip;

# Install SamTools
RUN cd /opt; wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2; \
    tar -vxjf samtools-1.10.tar.bz2; rm samtools-1.10.tar.bz2; \
    cd samtools-1.10; \
    ./configure; make; make install; \
    mv samtools /opt/bin;

# Install SamTools + BCFTools + htslib
RUN cd /opt; wget https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2; \
    tar -vxjf htslib-1.10.2.tar.bz2; rm htslib-1.10.2.tar.bz2; \
    cd htslib-1.10.2; \
    ./configure; make; make install; \
    mv htslib /opt/bin;
RUN cd /opt; wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2; \
    tar -vxjf samtools-1.10.tar.bz2; rm samtools-1.10.tar.bz2; \
    cd samtools-1.10; \
    ./configure; make; make install; \
    mv samtools /opt/bin;
RUN cd /opt; wget https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2; \
    tar -vxjf bcftools-1.10.2.tar.bz2; rm bcftools-1.10.2.tar.bz2; \
    cd bcftools-1.10.2; \
    ./configure; make; make install; \
    mv bcftools /opt/bin;

# Install Bedtools
RUN cd /opt; wget https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools-2.29.2.tar.gz; \
    tar xvfz bedtools-2.29.2.tar.gz; rm bedtools-2.29.2.tar.gz; \
    cd bedtools2; \
    make; \
    mv bin/* /opt/bin;

# Install FastQC
RUN cd /opt/bin; wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip; \
    unzip fastqc_v0.11.5.zip; rm fastqc_v0.11.5.zip;

# Install Picard
RUN cd /opt/bin; wget https://github.com/broadinstitute/picard/releases/download/2.21.7/picard.jar; \
    mv /opt/bin/* /usr/local/bin/; \
    chmod +x /usr/local/bin/picard.jar; \
    chmod 755 /usr/local/bin/FastQC/fastqc;
ENV PICARD "/usr/local/bin"
ENV PATH "$PATH:/usr/local/bin/FastQC"

# Install NGSCheckMate
RUN cd /opt/bin; wget https://github.com/parklab/NGSCheckMate/archive/v1.0.0.tar.gz; \
    tar xvfz v1.0.0.tar.gz; rm v1.0.0.tar.gz; \
    cd NGSCheckMate-1.0.0; \
    rm README.md LICENSE Documentation.pdf test_datafiles.txt; \
    mv * /usr/local/bin/;
ENV NCM_HOME="/usr/local/bin"

# Install MSIsensor-pro
RUN cd /opt/bin; git clone --branch 1.0.a https://github.com/xjtu-omics/msisensor-pro; \
    cd msisensor-pro/cpp; \
    make; \
    mv msisensor-pro /usr/local/bin/;

# Install BAM-ReadCount
RUN cd /opt; wget https://github.com/genome/bam-readcount/archive/v0.8.0.tar.gz; \
    tar xvfz v0.8.0.tar.gz; rm v0.8.0.tar.gz; \
    cd bam-readcount-0.8.0; \
    cmake .; make; \
    mv bin/bam-readcount /usr/local/bin/;

# Install Scripts
ENV repodir /seqprg
RUN mkdir -p /seqprg; \
    cd /seqprg; wget https://github.com/bcantarel/process_scripts/archive/version_1.1.5.tar.gz; tar xvfz version_1.1.5.tar.gz; \
    #git clone -b version_1.1.5 --single-branch --recursive https://github.com/bcantarel/process_scripts.git; \
    mv process_scripts-version_1.1.5 process_scripts; cd process_scripts; \
    cp -r alignment/* /usr/local/bin/; \
    cp -r variants/* /usr/local/bin/

ENV isdocker 1
ENV PATH "$PATH:/usr/local/bin"
ENV gitv "version_1.1.5"

WORKDIR /data/
