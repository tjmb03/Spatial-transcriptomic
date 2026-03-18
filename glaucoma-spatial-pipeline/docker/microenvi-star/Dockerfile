FROM ubuntu:22.04
LABEL maintainer="tjmb03"
LABEL version="1.0.0"
ARG STAR_VERSION=2.7.11a
ARG SAMTOOLS_VERSION=1.18
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
    wget curl gzip tar unzip \
    libncurses5-dev libbz2-dev liblzma-dev \
    zlib1g-dev gcc make python3 python3-pip \
    && rm -rf /var/lib/apt/lists/*
RUN wget -q https://github.com/alexdobin/STAR/releases/download/${STAR_VERSION}/STAR_${STAR_VERSION}.zip \
    && unzip -q STAR_${STAR_VERSION}.zip \
    && mv STAR_${STAR_VERSION}/Linux_x86_64_static/STAR /usr/local/bin/ \
    && chmod +x /usr/local/bin/STAR \
    && rm -rf STAR_${STAR_VERSION}*
RUN wget -q https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && cd samtools-${SAMTOOLS_VERSION} \
    && ./configure --prefix=/usr/local \
    && make -j4 && make install \
    && cd .. && rm -rf samtools-${SAMTOOLS_VERSION}*
RUN pip3 install --no-cache-dir numpy==1.26.3 pandas==2.1.3 scipy==1.11.4
RUN mkdir -p /pipeline/scripts /pipeline/config
RUN STAR --version && samtools --version | head -1
WORKDIR /pipeline
