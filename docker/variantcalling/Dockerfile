FROM broadinstitute/gatk:4.1.4.0
LABEL Author="Jeremy Mathews" \
      Maintainer="genome.school@gmail.com"

# install additional requirements
RUN apt-get update; apt-get install -y build-essential checkinstall; apt-get upgrade; \
    apt-get install -y --allow-unauthenticated libreadline-gplv2-dev libncursesw5-dev libssl-dev libsqlite3-dev tk-dev libgdbm-dev libc6-dev libbz2-dev zlib1g-dev liblzma-dev libcurl4-gnutls-dev libncurses5-dev wget unzip git default-jre default-jdk r-base cpanminus pigz


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
RUN cd /usr/src; wget https://bootstrap.pypa.io/pip/2.7/get-pip.py; \
    python2.7 get-pip.py;
RUN ln -sf /usr/local/bin/python2.7 /opt/miniconda/envs/gatk/bin/python

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

# Install Platypus
ENV C_INCLUDE_PATH $C_INCLUDE_PATH:/usr/local/include
ENV LIBRARY_PATH $LIBRARY_PATH:/usr/local/lib
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/usr/local/lib
RUN python -m pip install cython
RUN cd /usr/local/lib/python2.7/site-packages; \
    mv Cython cython;
RUN cd /opt; wget https://github.com/andyrimmer/Platypus/raw/master/release/AllReleases/Platypus_0.8.1.tgz; \
    tar -zxvf Platypus_0.8.1.tgz; rm Platypus_0.8.1.tgz; \
    cd Platypus_0.8.1; \
    ./buildPlatypus.sh; \
    cp -R * /opt/bin

# Install Bedtools
RUN cd /opt; wget https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools-2.29.2.tar.gz; \
    tar xvfz bedtools-2.29.2.tar.gz; rm bedtools-2.29.2.tar.gz; \
    cd bedtools2; \
    make; \
    cp -r bin/* /opt/bin;

# Install VCFTools
RUN cd /opt; wget https://github.com/vcftools/vcftools/releases/download/v0.1.14/vcftools-0.1.14.tar.gz; \
    tar xvfz vcftools-0.1.14.tar.gz; rm vcftools-0.1.14.tar.gz; \
    cd vcftools-0.1.14; \
    ./configure; make; make install;

# Install Picard + snpEff
RUN cd /opt/bin; wget https://github.com/broadinstitute/picard/releases/download/2.21.7/picard.jar;
RUN cd /opt/bin; wget https://sourceforge.net/projects/snpeff/files/snpEff_v4_3t_core.zip; \
    unzip snpEff_v4_3t_core.zip; rm snpEff_v4_3t_core.zip; \
    cd snpEff; \
    wget https://sourceforge.net/projects/snpeff/files/databases/v4_3/snpEff_v4_3_GRCh38.86.zip; \
    unzip snpEff_v4_3_GRCh38.86.zip; rm snpEff_v4_3_GRCh38.86.zip; \
    echo 'GRCh38.86.genome : Homo_sapiens' >> snpEff.config; \
    echo 'GRCh38.86.reference : ftp://ftp.ensembl.org/pub/release-92/gtf/' >> snpEff.config;
RUN cp -R /opt/bin/* /usr/local/bin/; \
    chmod +x /usr/local/bin/picard.jar; \
    chmod +x /usr/local/bin/snpEff/snpEff.jar; \
    chmod +x /usr/local/bin/snpEff/SnpSift.jar;
ENV PICARD "/usr/local/bin"
ENV SNPEFF_HOME "/usr/local/bin/snpEff"

# Install Strelka
RUN cd /opt; wget https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2; \
    tar -vxjf strelka-2.9.10.centos6_x86_64.tar.bz2; rm strelka-2.9.10.centos6_x86_64.tar.bz2; \
    mv strelka-2.9.10.centos6_x86_64 /usr/local/bin;

ENV PATH "$PATH:/usr/local/bin/strelka-2.9.10.centos6_x86_64/bin"

# Install Manta
RUN cd /opt; wget https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.release_src.tar.bz2; \
    tar -vxjf manta-1.6.0.release_src.tar.bz2; rm  manta-1.6.0.release_src.tar.bz2; \
    cd /usr/local/bin; \
    /opt/manta-1.6.0.release_src/configure --jobs=4; make -j4 install;

# Install Freebayes
RUN cd /opt; git clone --recursive --branch v1.2.0 git://github.com/ekg/freebayes.git; \
    cd freebayes; \
    make; make install;

# Install Parallel
RUN cd /opt/bin; wget http://ftp.gnu.org/gnu/parallel/parallel-latest.tar.bz2; \
    tar -vxjf parallel-latest.tar.bz2; rm parallel-latest.tar.bz2; \
    cd parallel*; \
    ./configure --prefix=/usr/local && make && make install;

# Install Vt
RUN cd /opt/bin; wget https://github.com/atks/vt/archive/0.57721.tar.gz; \
    tar xvfz 0.57721.tar.gz; rm -rf 0.57721.tar.gz; \
    cd vt-0.57721; \
    make; \
    mv vt /usr/local/bin/;

# Install Shimmer
RUN R -e "install.packages('statmod')"
RUN cpanm CPAN::Meta Module::Build
RUN cd /opt; wget https://github.com/bcantarel/Shimmer/archive/v0.1.2.zip; \
    unzip v0.1.2.zip; rm v0.1.2.zip; \
    cd  Shimmer-0.1.2; \
    perl Build.PL; \
    ./Build; ./Build test; ./Build install;

# Install Scripts
ENV repodir /seqprg
RUN mkdir -p /seqprg; \
    cd /seqprg; wget https://github.com/bcantarel/process_scripts/archive/version_1.1.5.tar.gz; tar xvfz version_1.1.5.tar.gz; \
    #git clone -b version_1.1.5 --single-branch --recursive https://github.com/bcantarel/process_scripts.git; \
    mv process_scripts-version_1.1.5 process_scripts; cd process_scripts; \
    cp -r alignment/* /usr/local/bin/; \
    cp -r variants/* /usr/local/bin/

RUN ln -sf /usr/local/bin/python2.7 /opt/miniconda/envs/gatk/bin/python
ENV isdocker 1
ENV PATH "$PATH:/usr/local/bin"

WORKDIR /data/
