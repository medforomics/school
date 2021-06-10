FROM r-base:4.0.2
LABEL Author="Jeremy Mathews" \
      Maintainer="genome.school@gmail.com"

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    libfftw3-dev libfftw3-doc \
    gcc \
    libcurl4-openssl-dev libssl-dev\
    libxml2-dev \
    git \
    procps \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install R Packages
RUN R -e 'install.packages(c("RColorBrewer","tidyverse","ggfortify","heatmaply","gplots","BiocManager"), dependencies=TRUE)'
RUN Rscript -e 'BiocManager::install(c("ballgown","qusage","edgeR", "DESeq2"), update=FALSE, quietly = TRUE)'

# Install Scripts
ENV repodir /seqprg
RUN mkdir -p /seqprg; \
    cd /seqprg; \
    git clone -b version_1.1.5 --single-branch --recursive https://github.com/bcantarel/process_scripts.git; \
    cd process_scripts; \
    cp -r genect_rnaseq/* /usr/local/bin/;

ENV isdocker 1
ENV PATH "$PATH:/usr/local/bin"

CMD []

WORKDIR /data/

