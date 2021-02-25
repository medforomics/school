FROM ubuntu:16.04
LABEL Author="Jeremy Mathews" \
      Maintainer="b.cantarel@gmail.com"

# install additional requirements
RUN apt-get update; apt-get install -y build-essential checkinstall; \
    apt-get install -y libreadline-gplv2-dev cmake patch libncursesw5-dev libssl-dev libsqlite3-dev tk-dev libgdbm-dev libc6-dev libbz2-dev zlib1g-dev liblzma-dev libcurl4-gnutls-dev libncurses5-dev wget unzip git default-jre default-jdk libtbb-dev pigz;


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
RUN python -m ensurepip --default-pip; \
    python -m pip install pysam argparse

# Install SamTools
RUN cd /opt; wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2; \
    tar -vxjf samtools-1.10.tar.bz2; rm samtools-1.10.tar.bz2; \
    cd samtools-1.10; \
    ./configure; make; make install; \
    cp samtools /opt/bin;

# Install HISAT2
RUN cd /opt; wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-source.zip; \
    unzip hisat2-2.1.0-source.zip; rm hisat2-2.1.0-source.zip; \
    cd hisat2-2.1.0; \
    make; \
    cp -r hisat* /usr/local/bin;

# Install Scripts
ENV repodir /seqprg
RUN mkdir -p /seqprg; \
    cd /seqprg; wget https://github.com/bcantarel/process_scripts/archive/version_1.0.9.tar.gz; tar xvfz version_1.0.9.tar.gz; \
    #git clone -b version_1.0.9 --single-branch --recursive https://github.com/bcantarel/process_scripts.git; \
    mv process_scripts-version_1.0.9 process_scripts; cd process_scripts; \
    cp -r alignment/* /usr/local/bin/; \
    cp -r variants/* /usr/local/bin/

ENV isdocker 1
ENV PATH "$PATH:/usr/local/bin"

WORKDIR /data/
