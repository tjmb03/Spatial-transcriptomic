FROM bioconductor/bioconductor_docker:RELEASE_3_18
LABEL maintainer="tjmb03"
LABEL version="1.0.0"

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    patch libgsl-dev liblapack-dev libblas-dev gfortran \
    libgomp1 libomp-dev g++ gcc make \
    && rm -rf /var/lib/apt/lists/*

RUN sed -i 's/^CXX =.*/CXX = g++ -std=c++17/' $(R RHOME)/etc/Makeconf && \
    sed -i 's/^CXXFLAGS =.*/CXXFLAGS = -O2 -fPIC -std=c++17/' $(R RHOME)/etc/Makeconf && \
    grep '^CXX ' $(R RHOME)/etc/Makeconf

COPY glaucoma-spatial-pipeline/docker/glaucoma-seurat/install.R /tmp/install.R
RUN Rscript /tmp/install.R

RUN Rscript -e "install.packages('BiocManager', repos='https://cloud.r-project.org')"

RUN Rscript -e "install.packages('BH', repos='https://cloud.r-project.org')"

RUN Rscript -e "BiocManager::install('spdep', ask=FALSE, update=FALSE, force=TRUE); stopifnot(requireNamespace('spdep', quietly=TRUE))"

RUN Rscript -e "BiocManager::install('org.Mm.eg.db', ask=FALSE, update=FALSE, force=TRUE); stopifnot(requireNamespace('org.Mm.eg.db', quietly=TRUE))"

RUN Rscript -e "BiocManager::install('DESeq2', ask=FALSE, update=FALSE, force=TRUE); stopifnot(requireNamespace('DESeq2', quietly=TRUE))"

RUN Rscript -e "BiocManager::install(c('limma','BiocParallel'), ask=FALSE, update=FALSE, force=TRUE)"

RUN mkdir -p /pipeline/scripts /pipeline/config
COPY glaucoma-spatial-pipeline/config/preregistration.py /pipeline/config/

WORKDIR /pipeline
CMD ["Rscript", "--help"]
