FROM ubuntu:16.04
LABEL Author="Jeremy Mathews" \
      Maintainer="b.cantarel@gmail.com"

# install additional requirements
RUN apt-get update; apt-get install -y build-essential checkinstall; \
    apt-get install -y libreadline-gplv2-dev libncursesw5-dev libssl-dev libsqlite3-dev tk-dev libgdbm-dev libc6-dev libbz2-dev zlib1g-dev liblzma-dev libcurl4-gnutls-dev libncurses5-dev wget unzip git default-jre default-jdk pigz;

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
    python -m pip install pysam

# Install BWA + BWAkit
RUN cd /usr/src; wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2; \
    tar -vxjf bwa-0.7.17.tar.bz2; rm bwa-0.7.17.tar.bz2; \
    cd bwa-0.7.17; make;

ENV PATH "$PATH:/usr/src/bwa-0.7.17"
RUN cd /usr/src; wget https://sourceforge.net/projects/bio-bwa/files/bwakit/bwakit-0.7.15_x64-linux.tar.bz2/download; \
    tar -vxjf download; rm download; \
    cd bwa.kit; \
    cp -R * /usr/local/bin;

# Install SamTools
RUN cd /opt; wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2; \
    tar -vxjf samtools-1.10.tar.bz2; rm samtools-1.10.tar.bz2; \
    cd samtools-1.10; \
    ./configure; make; make install; \
    cp samtools /opt/bin;

# Install Picard + FGBio
RUN cd /opt/bin; wget https://github.com/broadinstitute/picard/releases/download/2.21.7/picard.jar; \
    cd /opt/bin; wget https://github.com/fulcrumgenomics/fgbio/releases/download/1.1.0/fgbio-1.1.0.jar; \
    cp -R /opt/bin/* /usr/local/bin/; \
    chmod +x /usr/local/bin/picard.jar; \
    chmod +x /usr/local/bin/fgbio-1.1.0.jar;
ENV PICARD "/usr/local/bin"
RUN echo '#!/bin/bash' > /usr/local/bin/fgbio; \
    echo 'java -jar /usr/local/bin/fgbio-1.1.0.jar "$@"' >> /usr/local/bin/fgbio; \
    chmod +x /usr/local/bin/fgbio;

# Install Scripts
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
