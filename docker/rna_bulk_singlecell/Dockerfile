FROM quay.io/cellgeni/jupyter:r-base-20210222

LABEL Author="Brandi Cantarel" \
      Maintainer="genome.school@gmail.com"

# general OS packages
RUN apt-get update -y -qq 

# clean up
RUN apt-get clean && apt-get install -y libfftw3-dev libfftw3-doc &&\
    rm -rf /var/lib/apt/lists/*

## R packages
RUN Rscript -e 'install.packages(c("tidyverse", "rJava", "umap", "ggplot2", "ggfortify", "igraph","devtools", "Rmagic", "lsa", "uwot", "optparse", "Seurat","Ecfun","gridExtra","viridis","NMI","clustree","heatmaply"), dependencies = TRUE)'

RUN Rscript -e 'BiocManager::install(c("SummarizedExperiment", "SingleCellExperiment", "DropletUtils", "LoomExperiment", "Rhdf5lib", "scater", "scran", "RUVSeq", "sva", "MultiAssayExperiment", "batchelor", "BiocParallel", "SC3", "destiny", "pcaMethods", "zinbwave", "GenomicAlignments", "M3Drop", "switchde", "biomaRt", "Matrix.utils","SingleR","genefu","celldex","glmGamPoi"), update=FALSE, quietly = TRUE)'

RUN Rscript -e 'devtools::install_github(c("cellgeni/sceasy", "mojaveazure/loomR", "immunogenomics/harmony", "cole-trapnell-lab/leidenbase", "cole-trapnell-lab/monocle3","chris-mcginnis-ucsf/DoubletFinder"));'

RUN Rscript -e 'BiocManager::install(c("ballgown","qusage","edgeR", "DESeq2"), update=FALSE, quietly = TRUE)'

# clean conda cache
RUN  mamba clean --index-cache --tarballs --yes

# fix permissions
RUN fix-permissions /usr/lib/R/ && \
    fix-permissions /usr/local/lib/R/site-library

RUN echo "$image_name=$tag_name" >> /sanger/image.info && \ 
    conda list > /sanger/conda.info && \
    Rscript -e "as.data.frame(installed.packages()[,c(1,3:4)])[,1:2,drop=FALSE]" > /sanger/r-packages.info

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