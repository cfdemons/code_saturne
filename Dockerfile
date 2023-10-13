FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y && \
    apt-get install -y build-essential autoconf python3 gfortran libtool rsync openmpi-bin openmpi-common libopenmpi-dev libopenmpi3 libhdf5-openmpi-103 libhdf5-openmpi-dev liblapack-dev less pyqt5-dev python3-pyqt5 pyqt5-dev-tools&& \
    rm -rf /var/lib/apt/lists/*

RUN apt-get update -y && \
    apt-get install -y curl && \
    rm -rf /var/lib/apt/lists/*

RUN curl --create-dirs -O --output-dir /tmp/code_saturne \
    https://www.code-saturne.org/releases/code_saturne-7.0.2.tar.gz

RUN cd /tmp/code_saturne && tar xfz code_saturne-7.0.2.tar.gz

RUN mkdir /tmp/code_saturne/build

ARG MAKEFLAGS=--jobs=8

ARG INSTALL_DIR=/opt/code_saturne-7.0.2

WORKDIR /tmp/code_saturne/build

RUN ../code_saturne-7.0.2/configure \
    --prefix=$INSTALL_DIR \
    --disable-gui \
    --without-hdf5 \
    --without-cgns \
    --without-med \
    --without-metis \
    --without-scotch \
    --disable-static \
    --enable-debug \
    CXX=mpicxx \
    CC=mpicc \
    FC=mpif90 \
    CFLAGS="-g" \
    LDFLAGS="-g"

RUN make -j 8

RUN make install

ENV PATH=/opt/code_saturne-7.0.2/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

WORKDIR /opt

RUN mkdir /opt/packages

RUN tar cfz /opt/packages/code_saturne-7.0.2-jammy.tar.gz code_saturne-7.0.2
