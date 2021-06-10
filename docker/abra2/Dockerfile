FROM ubuntu:16.04
LABEL Author="Jeremy Mathews" \
      Maintainer="genome.school@gmail.com"

# install additional requirements
RUN apt-get update; apt-get install -y build-essential checkinstall; \
    apt-get install -y libreadline-gplv2-dev cmake patch libncursesw5-dev libssl-dev libsqlite3-dev tk-dev libgdbm-dev libc6-dev libbz2-dev zlib1g-dev liblzma-dev libcurl4-gnutls-dev libncurses5-dev wget unzip git default-jre default-jdk libtbb-dev pigz;

ARG ABRA2_VERSION=2.20
ENV ABRA2_VERSION ${ABRA2_VERSION}
ENV JAVA_OPTS "-Xmx16G"

RUN mkdir -p /opt;
RUN mkdir -p /opt/bin;

RUN cd /opt/bin; wget https://github.com/mozack/abra2/releases/download/v${ABRA2_VERSION}/abra2-${ABRA2_VERSION}.jar; chmod 755 /opt/bin/abra2-${ABRA2_VERSION}.jar; ln -s /opt/bin/abra2-${ABRA2_VERSION}.jar /usr/local/bin/abra2.jar 

# Install SamTools
RUN cd /opt; wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2; \
    tar -vxjf samtools-1.10.tar.bz2; rm samtools-1.10.tar.bz2; \
    cd samtools-1.10; \
    ./configure; make; make install; \
    cp samtools /usr/local/bin;

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

WORKDIR /data/
