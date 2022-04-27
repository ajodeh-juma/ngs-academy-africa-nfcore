# base image
FROM ubuntu:xenial

# metadata
LABEL base.image="ubuntu:xenial"
LABEL authors="@ajodeh-juma" 
LABEL description="Docker image containing all software requirements for the simple RNA-seq pipeline"


# install software-properties-common for add-apt-repository
RUN apt-get update && apt-get install -y \
    software-properties-common


# add universe repository for texlive-xetex and texlive-math-extra
RUN add-apt-repository universe

# install dependencies; cleanup apt garbage
RUN apt-get update && apt-get install -y \
    unzip \
    wget \
    curl \
    perl \
    default-jre \
    build-essential \
    zlib1g-dev \
    libncurses5-dev libgdbm-dev libnss3-dev libssl-dev libreadline-dev \
    libffi-dev libbz2-dev


# install python3
RUN wget https://www.python.org/ftp/python/3.9.10/Python-3.9.10.tar.xz && \
    tar xf Python-3.9.10.tar.xz && \
    cd Python-3.9.10 && \
    ./configure --enable-optimizations --with-lto --enable-ipv6 --with-tzpath=/usr/share/zoneinfo && \
    make install && \
    make maninstall

# install multiqc
RUN pip3 install "multiqc==1.12"

#  # install fastqc
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
    unzip fastqc_v0.11.9.zip && \
    rm fastqc_v0.11.9.zip && \
    chmod +x FastQC/fastqc && \
    mkdir /data

# set PATH 
ENV PATH="${PATH}:/FastQC/"


# install salmon
RUN curl -sSL https://github.com/COMBINE-lab/salmon/releases/download/v1.8.0/salmon-1.8.0_linux_x86_64.tar.gz | tar xz \
 && mv /salmon-*/bin/* /usr/bin/ \
 && mv /salmon-*/lib/* /usr/lib/

# set working directory
WORKDIR /data