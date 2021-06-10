FROM ubuntu:16.04
LABEL Author="Jeremy Mathews" \
      Maintainer="genome.school@gmail.com"


#install additional requirements
RUN apt-get update; apt-get install -y build-essential checkinstall; \
    apt-get install -y libreadline-gplv2-dev libncursesw5-dev libssl-dev libsqlite3-dev tk-dev libgdbm-dev libc6-dev libbz2-dev zlib1g-dev liblzma-dev libcurl4-gnutls-dev libncurses5-dev libffi-dev wget unzip git default-jre default-jdk cmake gfortran hdf5-tools libboost-date-time-dev libboost-program-options-dev libboost-system-dev libboost-filesystem-dev libboost-iostreams-dev libhdf5-dev fort77 xorg-dev libblas-dev gcc-multilib gobjc++ aptitude texlive-latex-base libcairo2-dev pigz;
    
RUN cd /usr/src; wget https://www.python.org/ftp/python/3.8.1/Python-3.8.1.tgz; \
    tar -zxvf Python-3.8.1.tgz; rm Python-3.8.1.tgz; \
    cd Python-3.8.1; \
    ./configure --enable-optimizations; \
    make altinstall;
RUN cd /usr/local/bin; \
    mv python3.8 python;
RUN python -m ensurepip --default-pip;

# Make destination and work directory
RUN mkdir -p /opt;
RUN mkdir -p /opt/bin;
#ENV http_proxy http://proxy.swmed.edu:3128/
#ENV https_proxy http://proxy.swmed.edu:3128/

# Install Subread
RUN cd /opt; wget https://sourceforge.net/projects/subread/files/subread-2.0.0/subread-2.0.0-Linux-x86_64.tar.gz; \
    tar zxvf subread-2.0.0-Linux-x86_64.tar.gz; rm subread-2.0.0-Linux-x86_64.tar.gz; \
    cd subread-2.0.0-Linux-x86_64/bin; \
    cp -r * /usr/local/bin;

# Install StringTie
RUN cd /opt; wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.1.1.Linux_x86_64.tar.gz; \
    tar xvfz stringtie-2.1.1.Linux_x86_64.tar.gz; rm stringtie-2.1.1.Linux_x86_64.tar.gz; \
    cd stringtie-2.1.1.Linux_x86_64; \
    cp * /usr/local/bin;

RUN cd /opt; wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2; \
    tar -vxjf samtools-1.10.tar.bz2; rm samtools-1.10.tar.bz2; \
    cd samtools-1.10; \
    ./configure; make; make install; \
    cp samtools /usr/local/bin/;

RUN cd /opt; wget https://github.com/griffithlab/regtools/archive/0.5.2.tar.gz; tar xvfz 0.5.2.tar.gz; cd regtools-0.5.2/; mkdir build; cd build/; cmake ..; make; cp regtools /usr/local/bin/

# Install Scripts
RUN cd /opt; wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz; \
    tar xvfz bedtools-2.30.0.tar.gz; rm bedtools-2.30.0.tar.gz; \
    cd bedtools2; \
    make; \
    cp bin/* /usr/local/bin;

ENV repodir /seqprg
RUN mkdir -p /seqprg; \
    cd /seqprg; wget https://github.com/bcantarel/process_scripts/archive/version_1.1.5.tar.gz; tar xvfz version_1.1.5.tar.gz; \
    #git clone -b version_1.1.5 --single-branch --recursive https://github.com/bcantarel/process_scripts.git; \
    mv process_scripts-version_1.1.5 process_scripts; cd process_scripts; \
    cp -r genect_rnaseq/* /usr/local/bin/

ENV isdocker 1
ENV PATH "$PATH:/usr/local/bin"
WORKDIR /data/
