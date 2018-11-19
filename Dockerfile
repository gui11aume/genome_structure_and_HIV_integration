FROM ubuntu:14.04.5

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
            r-base          \
            r-recommended   \
            tabix           \
            pigz            \
            cmake           \
    && apt-get clean &&     \
    rm -rf /var/lib/apt/lists*

# bwa
RUN cd / && \
    git clone https://github.com/lh3/bwa && \
    cd bwa && \
    git checkout 5961611c358e480110793bbf241523a3cfac049b && \
    make && \
    cp bwa /usr/local/bin

# zerone
RUN cd / && \
    git clone https://github.com/nanakiksc/zerone && \
    cd zerone && \
    git checkout 449e8289244b38185f0ff37472d3ff55288d9615 && \
    make && \
    cp zerone /usr/local/bin

# starcode
RUN cd / && \
    git clone https://github.com/gui11aume/starcode && \
    cd starcode && \
    git checkout 2255df95d51536be70708c65dc7c5f83bda4be7a && \
    make && \
    ln -s /starcode/starcode /usr/local/bin/starcode

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

# samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 && \
    tar -xvf samtools-1.3.1.tar.bz2 && \
    cd samtools-1.3.1 && \
    ./configure && \
    make && make install

# cooler (from pip)
RUN pip install cooler

# R libraries
RUN R -e 'source("http://bioconductor.org/biocLite.R"); biocLite("GenomicRanges")'