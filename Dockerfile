FROM ubuntu:14.04.5

# Add repositories
RUN apt-get update -qq &&   \
    apt-get install -y      \
            apt-transport-https

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
    echo "deb https://cloud.r-project.org/bin/linux/ubuntu trusty-cran35/" >> /etc/apt/sources.list

# Software from trusty repositories
RUN apt-get update -qq &&   \
    apt-get install -y      \
            build-essential \
            git             \
            wget            \
            gzip            \
            python          \
            python-dev      \
            python-pip      \
            python-numpy    \
            python-scipy    \
            python-pandas   \
            python-h5py     \
            zlib1g-dev      \
            libxml2-dev     \
            libmagic-dev    \
            libhdf5-dev     \
            libbz2-dev      \
            liblzma-dev     \
            libncurses5-dev \
            libfuse-dev     \
	    libcurl4-openssl-dev \
	    libxml2-dev     \
            r-base-core     \
            r-base-dev      \
	    r-cran-mgcv	    \
            tabix           \
            pigz            \
            cmake           \
    && apt-get clean &&     \
    rm -rf /var/lib/apt/lists/*

# bwa
RUN cd / && \
    git clone https://github.com/lh3/bwa && \
    cd bwa && \
    make && \
    cp bwa /usr/local/bin

# zerone
RUN cd / && \
    git clone https://github.com/nanakiksc/zerone && \
    cd zerone && \
    make && \
    cp zerone /usr/local/bin

# starcode
RUN cd / && \
    git clone https://github.com/gui11aume/starcode && \
    cd starcode && \
    make && \
    ln -s /starcode/starcode /usr/local/bin/starcode

# pandoc (for rmarkdown)
RUN cd / && \
    wget https://github.com/jgm/pandoc/releases/download/2.5/pandoc-2.5-1-amd64.deb && \
    dpkg -i pandoc-2.5-1-amd64.deb && \
    rm pandoc-2.5-1-amd64.deb

# ncbi tools
RUN cd / && \
    mkdir ncbi && \
    cd ncbi && \
    git clone http://github.com/ncbi/ngs && \
    git clone http://github.com/ncbi/ncbi-vdb && \
    git clone http://github.com/ncbi/sra-tools && \
    ./ngs/configure && make -C ngs/ngs-sdk install && \
    ./ncbi-vdb/configure && make -C ncbi-vdb install && \
    ./sra-tools/configure && make -C sra-tools install && \
    mv /usr/local/ncbi/sra-tools/bin/* /usr/local/bin/

# # samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 && \
    tar -xvf samtools-1.3.1.tar.bz2 && \
    cd samtools-1.3.1 && \
    ./configure && \
    make && make install

# cooler (from pip)
RUN pip install cooler

# sklearn (from pip)
RUN pip install sklearn

# R libraries
RUN R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install(c("GenomicRanges","Sushi","data.table","digest","stringr","pheatmap"), version = "3.8");'
RUN R -e 'install.packages("gplots", repo="https://cran.cnr.berkeley.edu/")'
RUN R -e 'install.packages("ggplot2", repo="https://cran.cnr.berkeley.edu/")'
RUN R -e 'install.packages("RColorBrewer", repo="https://cran.cnr.berkeley.edu/")'
RUN R -e 'install.packages("gridExtra", repo="https://cran.cnr.berkeley.edu/")'
RUN R -e 'install.packages("scales", repo="https://cran.cnr.berkeley.edu/")'
RUN R -e 'install.packages("rmarkdown", repo="https://cran.cnr.berkeley.edu/")'
RUN R -e 'install.packages("tidyr", repo="https://cran.cnr.berkeley.edu/")'