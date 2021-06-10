FROM etal/cnvkit:0.9.7
LABEL Author="Jeremy Mathews" \
      Maintainer="genome.school@gmail.com"

USER root
# install additional requirements
RUN apt-get update; apt-get install -y build-essential checkinstall; \
    apt-get install -y libreadline-gplv2-dev libncursesw5-dev libssl-dev libsqlite3-dev tk-dev libgdbm-dev libc6-dev libbz2-dev zlib1g-dev liblzma-dev libcurl4-gnutls-dev libncurses5-dev libffi-dev wget unzip git default-jre default-jdk cmake gfortran hdf5-tools libboost-date-time-dev libboost-program-options-dev libboost-system-dev libboost-filesystem-dev libboost-iostreams-dev libhdf5-dev fort77 xorg-dev libblas-dev gcc-multilib gobjc++ aptitude texlive-latex-base libcairo2-dev pigz;

# Make destination and work directory
RUN mkdir -p /opt;
RUN mkdir -p /opt/bin;
#ENV http_proxy http://proxy.swmed.edu:3128/
#ENV https_proxy http://proxy.swmed.edu:3128/

# Install Python
RUN cd /usr/src; wget https://www.python.org/ftp/python/3.8.1/Python-3.8.1.tgz; \
    tar -zxvf Python-3.8.1.tgz; rm Python-3.8.1.tgz; \
    cd Python-3.8.1; \
    ./configure --enable-optimizations; \
    make altinstall;
RUN cd /usr/local/bin; \
    mv python3.8 python;
RUN python -m ensurepip --default-pip;

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
    cp bin/* /opt/bin;

# Install VCFTools
RUN cd /opt; wget https://github.com/vcftools/vcftools/releases/download/v0.1.14/vcftools-0.1.14.tar.gz; \
    tar xvfz vcftools-0.1.14.tar.gz; rm vcftools-0.1.14.tar.gz; \
    cd vcftools-0.1.14; \
    ./configure; make; make install;

# Install Picard + snpEff
RUN cd /opt/bin; wget https://sourceforge.net/projects/snpeff/files/snpEff_v4_3t_core.zip; \
    unzip snpEff_v4_3t_core.zip; rm snpEff_v4_3t_core.zip; \
    cd snpEff; \
    wget https://sourceforge.net/projects/snpeff/files/databases/v4_3/snpEff_v4_3_GRCh38.86.zip; \
    unzip snpEff_v4_3_GRCh38.86.zip; rm snpEff_v4_3_GRCh38.86.zip; \
    echo 'GRCh38.86.genome : Homo_sapiens' >> snpEff.config; \
    echo 'GRCh38.86.reference : ftp://ftp.ensembl.org/pub/release-92/gtf/' >> snpEff.config;
RUN cp -R /opt/bin/* /usr/local/bin/; \
    chmod +x /usr/local/bin/snpEff/snpEff.jar; \
    chmod +x /usr/local/bin/snpEff/SnpSift.jar;
ENV SNPEFF_HOME "/usr/local/bin/snpEff"

# Install Pindel
RUN git clone https://github.com/xjtu-omics/pindel; \
    cd pindel; \
    ./INSTALL /opt/samtools-1.10/htslib-1.10; \
    cp *pindel* /usr/local/bin/;
    
# Install Delly2
RUN cd /opt/bin; git clone --branch v0.8.2 --recursive https://github.com/dellytools/delly.git; \
    cd delly; \
    make all; \
    cp src/* /usr/local/bin/; \
    mv /usr/local/bin/delly /usr/local/bin/delly2;

# Install SvABA
RUN cd /opt/bin; git clone --branch 1.1.0 --recursive https://github.com/walaj/svaba; \
    cd svaba; \
    ./configure; \
    make && make install; \
    cp bin/* /usr/local/bin/;

# Install ITDSeek
RUN cd /opt/bin; wget https://github.com/bcantarel/itdseek/archive/v1.3.2.tar.gz; \
    tar xvfz v1.3.2.tar.gz; rm v1.3.2.tar.gz; \
    cd itdseek-1.3.2; \
    mv itd* /usr/local/bin/;


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
ENV SNPEFF_HOME "/usr/local/bin/snpEff"
WORKDIR /data/
