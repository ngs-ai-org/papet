# building in ubuntu 20.04
FROM ubuntu@sha256:b795f8e0caaaacad9859a9a38fe1c78154f8301fdaf0872eaf1520d66d9c0b98 AS builder
SHELL ["/bin/bash", "-c"]

ENV CONTAINER_NAME="NGS-AI PAPET"
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
# ENV PATH="/usr/local/bin:$PATH"

# turn off interactive program installation
ARG DEBIAN_FRONTEND=noninteractive

# version numbers
ARG boost_version="1_78_0"
ARG htslib_version="1.13" 
ARG samtools_version=${htslib_version}
ARG meson_version="1.1.0"
ARG ninja_version="1.11.1"
ARG pbcopper_version="v2.0.0"
ARG pbbam_version="v2.0.0"
ARG gtest_version="release-1.11.0"

# install general programs
RUN apt update && apt install -y wget && apt install -y 'python3' 'python3-pip' 'git' 'cmake' && \
# symbolic link to python
    ln -s /usr/bin/python3 /usr/bin/python && \
# install gcc and g++ v8 and make them by default
    apt update && apt install -y 'g++-8' 'gcc-8' && \
	ln -fs /usr/bin/gcc-8 /usr/bin/gcc && \
	ln -fs /usr/bin/g++-8 /usr/bin/g++ && \
## install building tools
    python3 -m pip install meson==${meson_version} && \
    python3 -m pip install ninja==${ninja_version} && \
    apt install -y autoconf automake make && \
## install dependencies
    ## HTSlib
    apt install -y perl && \
    ## HTSlib
    apt install -y zlib1g-dev && \
    ## HTSlib, boost
    apt install -y libbz2-dev && \
    ## HTSlib
    apt install -y liblzma-dev && \
    ## HTSlib
    apt install -y libcurl4-gnutls-dev && \
    ## HTSlib
    apt install -y libcurl4-gnutls-dev && \
    ## HTSlib
    apt install -y libncurses5-dev && \
    ## boost
    apt install -y build-essential && \
    ## boost
    apt install -y autotools-dev && \
    ## boost
    apt install -y libicu-dev && \
    ## pbcopper and pbbam
    apt install -y pkg-config

# install HTSlib
RUN wget https://github.com/samtools/htslib/releases/download/${htslib_version}/htslib-${htslib_version}.tar.bz2 -O htslib.tar.bz2  && \
    tar -xjvf htslib.tar.bz2  && \
    cd htslib-${htslib_version}  && \
## build
    make  && \
## install
    make install prefix=/usr && \
## clean
    cd .. && \
    rm -rf htslib*

# install samtools
RUN wget https://github.com/samtools/samtools/releases/download/${samtools_version}/samtools-${samtools_version}.tar.bz2 -O samtools.tar.bz2 && \
    tar -xjvf samtools.tar.bz2 && \
    cd samtools-${samtools_version} && \
## build
    ./configure --with-libhts=/usr/lib/libhts.so && \
    make && \
## install
    make prefix=/usr install && \
## clean
    cd .. && \
    rm -rf samtools*

# install boost
RUN wget -O boost_${boost_version}.tar.gz https://boostorg.jfrog.io/artifactory/main/release/1.78.0/source/boost_${boost_version}.tar.gz && \
    tar -zvxf boost_${boost_version}.tar.gz && \
## build
    cd boost_${boost_version}/ && ./bootstrap.sh --prefix=/usr && \
## install
    ./b2 install && \
## clean
   cd .. && \
   rm -rf boost_${boost_version} boost_${boost_version}.tar.gz

# install GoogleTest
## clone
RUN git clone --branch ${gtest_version} https://github.com/google/googletest.git && \
## compile and install in /usr/lib
	cd googletest/ && \
	cmake -DCMAKE_INSTALL_RPATH=/usr/lib -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE . && \
	make install && \
## clean
	cd .. && \
	rm -rf googletest/

# install pbcopper
RUN git clone --branch ${pbcopper_version} https://github.com/PacificBiosciences/pbcopper.git && \
    cd pbcopper/ && \
    mkdir -p build/ && \
    cd build/ && \
## build
    meson setup --prefix /usr .. && \
    ninja && \
## install
    ninja install && \
## clean
    cd ../.. && \
    rm -rf pbcopper*

# install pbbam
RUN git clone --branch ${pbbam_version} https://github.com/PacificBiosciences/pbbam.git && \
    cd pbbam/ && \
    mkdir -p build/ && \
    cd build/ && \
## build
    meson setup --prefix /usr .. && \
    ninja && \
## install
    ninja install && \

## clean
    cd ../.. && \
    rm -rf pbbam

# install ngsaipp
RUN git clone https://github.com/ngs-ai-org/ngsaipp.git && \
    cd ngsaipp/ && \
## build
    cmake . && \
    make && \
## install
    make install && \
## clean
    cd .. && \
    rm -rf ngsaipp*/

# install papet
RUN git clone https://github.com/ngs-ai-org/papet.git && \
    cd papet && \
## build
    cmake . && \
    make && \
## install
    make install && \
## clean
    cd .. && \
    rm -rf papet*/

# clean
RUN apt clean && apt-get clean

# final image
FROM ubuntu@sha256:b795f8e0caaaacad9859a9a38fe1c78154f8301fdaf0872eaf1520d66d9c0b98
SHELL ["/bin/bash", "-c"]
ENV CONTAINER_NAME="NGS-AI PAPET"
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
# turn off interactive program installation
ARG DEBIAN_FRONTEND=noninteractive
# htslib
COPY --from=builder /usr/lib/libhts* /usr/lib/
COPY --from=builder /usr/include/htslib /usr/include/
COPY --from=builder /lib/x86_64-linux-gnu/libcurl-gnutls.so.4 /lib/x86_64-linux-gnu/
# samtools
COPY --from=builder /usr/bin/*sam* /usr/bin/
# COPY --from=builder /lib/x86_64-linux-gnu/libcurl-gnutls.so.4 /lib/x86_64-linux-gnu/
# boost
COPY --from=builder /usr/lib/libboost_program_options.so.1.78.0 /usr/lib/
COPY --from=builder /usr/lib/libboost_serialization.so.1.78.0 /usr/lib/
COPY --from=builder /usr/include/boost /usr/include/
COPY --from=builder /lib/x86_64-linux-gnu/libicui18n.so.66 /lib/x86_64-linux-gnu/
COPY --from=builder /lib/x86_64-linux-gnu/libicuuc.so.66   /lib/x86_64-linux-gnu/
# gmock
COPY --from=builder /usr/local/lib/libgmock* /usr/local/lib/
COPY --from=builder /usr/local/include/gmock/ /usr/local/include/
# gtest
COPY --from=builder /usr/local/lib/libgtest* /usr/local/lib/
COPY --from=builder /usr/local/include/gtest/ /usr/local/include/
# pbcopper
COPY --from=builder /usr/lib/x86_64-linux-gnu/libpbcopper* /usr/lib/x86_64-linux-gnu/
COPY --from=builder /usr/include/pbcopper/ /usr/include/
# pbbam
COPY --from=builder /usr/lib/x86_64-linux-gnu/libpbbam* /usr/lib/x86_64-linux-gnu/
COPY --from=builder /usr/include/pbbam/ /usr/include/
COPY --from=builder /lib/x86_64-linux-gnu/libcurl-gnutls.so.4 /lib/x86_64-linux-gnu/
# ngsaipp
COPY --from=builder /usr/local/lib/libngsaipp* /usr/local/lib/
COPY --from=builder /usr/local/include/ngsaipp/ /usr/local/include/
# COPY --from=builder /lib/x86_64-linux-gnu/libcurl-gnutls.so.4 /lib/x86_64-linux-gnu/
# papet
COPY --from=builder /usr/local/bin/papet /usr/local/bin/
COPY --from=builder /lib/x86_64-linux-gnu/libnghttp2.so.14 /lib/x86_64-linux-gnu/
COPY --from=builder /lib/x86_64-linux-gnu/librtmp.so.1 /lib/x86_64-linux-gnu/
COPY --from=builder /lib/x86_64-linux-gnu/libssh.so.4 /lib/x86_64-linux-gnu/
COPY --from=builder /lib/x86_64-linux-gnu/libpsl.so.5 /lib/x86_64-linux-gnu/
COPY --from=builder /lib/x86_64-linux-gnu/libgssapi_krb5.so.2 /lib/x86_64-linux-gnu/
COPY --from=builder /lib/x86_64-linux-gnu/libldap_r-2.4.so.2 /lib/x86_64-linux-gnu/
COPY --from=builder /lib/x86_64-linux-gnu/liblber-2.4.so.2 /lib/x86_64-linux-gnu/
COPY --from=builder /lib/x86_64-linux-gnu/libbrotlidec.so.1 /lib/x86_64-linux-gnu/
COPY --from=builder /lib/x86_64-linux-gnu/libcrypto.so.1.1 /lib/x86_64-linux-gnu/
COPY --from=builder /lib/x86_64-linux-gnu/libkrb5.so.3 /lib/x86_64-linux-gnu/
COPY --from=builder /lib/x86_64-linux-gnu/libk5crypto.so.3 /lib/x86_64-linux-gnu/
COPY --from=builder /lib/x86_64-linux-gnu/libkrb5support.so.0 /lib/x86_64-linux-gnu/
COPY --from=builder /lib/x86_64-linux-gnu/libsasl2.so.2 /lib/x86_64-linux-gnu/
COPY --from=builder /lib/x86_64-linux-gnu/libgssapi.so.3 /lib/x86_64-linux-gnu/
COPY --from=builder /lib/x86_64-linux-gnu/libbrotlicommon.so.1 /lib/x86_64-linux-gnu/
COPY --from=builder /lib/x86_64-linux-gnu/libkeyutils.so.1 /lib/x86_64-linux-gnu/
COPY --from=builder /lib/x86_64-linux-gnu/libheimntlm.so.0 /lib/x86_64-linux-gnu/
COPY --from=builder /lib/x86_64-linux-gnu/libasn1.so.8 /lib/x86_64-linux-gnu/
COPY --from=builder /lib/x86_64-linux-gnu/libroken.so.18 /lib/x86_64-linux-gnu/
COPY --from=builder /lib/x86_64-linux-gnu/libkrb5.so.26 /lib/x86_64-linux-gnu/
COPY --from=builder /lib/x86_64-linux-gnu/libwind.so.0 /lib/x86_64-linux-gnu/
COPY --from=builder /lib/x86_64-linux-gnu/libhcrypto.so.4 /lib/x86_64-linux-gnu/
COPY --from=builder /lib/x86_64-linux-gnu/libheimbase.so.1 /lib/x86_64-linux-gnu/
COPY --from=builder /lib/x86_64-linux-gnu/libhx509.so.5 /lib/x86_64-linux-gnu/
COPY --from=builder /lib/x86_64-linux-gnu/libsqlite3.so.0 /lib/x86_64-linux-gnu/
