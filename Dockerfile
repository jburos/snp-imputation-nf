FROM ubuntu:16.04

# Update the repository sources list and install samtools package
RUN apt-get update
RUN apt-get -y install libcurl4-gnutls-dev
RUN apt-get -y install libcurl4-openssl-dev
RUN apt-get -y install libxml2-dev \
                        libssl-dev \
                        r-base \
                        r-base-dev \
                        default-jdk \
                        samtools \
                        git \
                        libboost-all-dev \
                        wget \
                        unzip \
                        htop \
                        sudo \
                        tabix \
                        curl \
                        build-essential \
                        python \
                        python3 \
                        python-pip \
                        python3-pip

RUN apt-get -y install bcftools vcftools

#install java
RUN apt-get -y install default-jre

#install nextflow
WORKDIR /usr/local/bin
RUN curl -o nextflow -fsSL get.nextflow.io
RUN chmod +x nextflow
RUN /usr/local/bin/nextflow

#install node
RUN curl -sL https://deb.nodesource.com/setup_7.x | sudo -E bash -
RUN apt-get install -y nodejs

#install awscli
RUN pip install awscli --upgrade

#install gsutil
RUN pip install gsutil

RUN mkdir /install

WORKDIR /install

# Get plink 1.9
#install plink
RUN wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20181202.zip && \
  unzip plink_linux_x86*.zip && \
  rm *.zip && \
  mv plink plink-1.9
ENV PATH /install:$PATH

#get plink 1.7 & any other packages we want from aptitude
RUN sudo apt-get update && \
  sudo DEBIAN_FRONTEND=noninteractive apt-get install -yq plink && \
  rm -rf /var/lib/apt/lists/*

#install shapeit
RUN wget https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r837.GLIBCv2.12.Linux.static.tgz && \
   tar -zxvf shapeit.v2.r837.GLIBCv2.12.Linux.static.tgz
ENV PATH /install/bin:$PATH

#install impute2
RUN wget https://mathgen.stats.ox.ac.uk/impute/impute_v2.3.2_x86_64_static.tgz && \
  tar -zxvf impute_v2.3.2_x86_64_static.tgz
ENV PATH /install/impute_v2.3.2_x86_64_static:$PATH

#Get gtools
RUN wget http://www.well.ox.ac.uk/~cfreeman/software/gwas/gtool_v0.7.5_x86_64.tgz && \
  tar zxvf gtool_v0.7.5_x86_64.tgz

RUN mkdir /workdir
WORKDIR /workdir

