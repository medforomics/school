FROM r-base:4.0.2
LABEL Author="Jeremy Mathews" \
      Maintainer="b.cantarel@gmail.com"

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    libfftw3-dev \
    gcc \
    libcurl4-openssl-dev \
    libxml2-dev \
    git \
    procps \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

#Install BIOCONDUCTOR Packages
RUN R -e 'install.packages("BiocManager")'
RUN R -e 'BiocManager::install("edgeR")'
RUN R -e 'BiocManager::install("DESeq2")'
RUN R -e 'BiocManager::install("qusage")'
RUN R -e 'BiocManager::install("ballgown")'

# Install R Packages
RUN R -e "install.packages('RColorBrewer', dependencies=TRUE, repos='http://cran.rstudio.com/')"; \
    R -e "install.packages('gplots', dependencies=TRUE, repos='http://cran.rstudio.com/')";

# Install Scripts
ENV repodir /seqprg
RUN mkdir -p /seqprg; \
    cd /seqprg; \
    git clone -b version_1.0.9 --single-branch --recursive https://github.com/bcantarel/process_scripts.git; \
    cd process_scripts; \
    cp -r genect_rnaseq/* /usr/local/bin/;

ENV isdocker 1
ENV PATH "$PATH:/usr/local/bin"

CMD []

WORKDIR /data/