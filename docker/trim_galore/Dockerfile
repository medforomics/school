FROM ubuntu:16.04
LABEL Author="Jeremy Mathews" \
      Maintainer="genome.school@gmail.com"

# install additional requirements
RUN apt-get update && apt-get install -y \
   wget \
   unzip \
   python3-pip \
   git

# Make destination and work directory
RUN mkdir -p /opt;
RUN mkdir -p /opt/bin;
#ENV http_proxy http://proxy.swmed.edu:3128/
#ENV https_proxy http://proxy.swmed.edu:3128/

# Install cutadapt
RUN pip3 install cutadapt==1.9.1

# Install trim_galore
RUN cd /opt; wget https://github.com/ComputationalSystemsBiology/EoulsanDockerFiles/raw/master/TrimAdapt/trim_galore_v0.4.1.zip; \
    unzip trim_galore_v0.4.1.zip; rm trim_galore_v0.4.1.zip; \
    ln -s /opt/trim_galore_zip/trim_galore /usr/local/bin/trim_galore;

# Install Scripts
ENV repodir /seqprg
RUN mkdir -p /seqprg; \
    cd /seqprg; wget https://github.com/bcantarel/process_scripts/archive/version_1.1.5.tar.gz; tar xvfz version_1.1.5.tar.gz; \
    #git clone -b version_1.1.5 --single-branch --recursive https://github.com/bcantarel/process_scripts.git; \
    mv process_scripts-version_1.1.5 process_scripts; cd process_scripts; \
    cp -r preproc_fastq/* /usr/local/bin/

ENV isdocker 1
ENV PATH "$PATH:/usr/local/bin"

WORKDIR /data/
