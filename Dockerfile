FROM fedora:43 AS build

# ── System packages ──────────────────────────────────────────────────────────
RUN dnf install -y \
        gcc gcc-c++ make cmake \
        gsl gsl-devel \
        hdf5 hdf5-devel \
        gdal gdal-devel \
        libgeotiff libgeotiff-devel \
        libtiff libtiff-devel \
        proj proj-devel \
        git wget unzip \
    && dnf clean all

# ── Clone external dependencies ──────────────────────────────────────────────
RUN git clone https://bitbucket.org/StevenHancock/tools.git /opt/tools

RUN mkdir -p /opt/cmpfit && \
    cd /tmp && \
    wget -q https://www.physics.wisc.edu/~craigm/idl/down/cmpfit-1.2.tar.gz && \
    tar -xzf cmpfit-1.2.tar.gz && \
    cp -a cmpfit-1.2/* /opt/cmpfit/ && \
    rm -rf /tmp/cmpfit-1.2*

# ── Extract bundled Riegl libraries ──────────────────────────────────────────
COPY riegl_libs /tmp/riegl_libs

RUN mkdir -p /opt/rivlib /opt/rdblib && \
    cd /tmp/riegl_libs && \
    tar -xzf rdblib-2.4.1-x86_64-linux.tar.gz -C /opt/rdblib --strip-components=1 && \
    unzip -qo rivlib-2_6_0-x86_64-linux-gcc11.zip -d /opt/rivlib && \
    rm -rf /tmp/riegl_libs

RUN echo "/opt/rivlib/rivlib-2_6_0-x86_64-linux-gcc11/lib" > /etc/ld.so.conf.d/rivlib.conf && \
    echo "/opt/rdblib/library"                             >> /etc/ld.so.conf.d/rivlib.conf && \
    ldconfig

# ── Copy source and build ───────────────────────────────────────────────────
COPY libclidar    /src/libclidar
COPY gedisimulator /src/gedisimulator
COPY voxelate     /src/voxelate
COPY voxel_lidar  /src/voxel_lidar
COPY salcamake    /src/salcamake
COPY rxptobinary  /src/rxptobinary

ENV ARCH=x86_64 \
    HANCOCKTOOLS_ROOT=/opt/tools \
    CMPFIT_ROOT=/opt/cmpfit \
    LIBCLIDAR_ROOT=/src/libclidar \
    VOXELATE_ROOT=/src/voxelate \
    GEDIRAT_ROOT=/src/gedisimulator \
    GSL_ROOT=/usr/lib64 \
    HDF5_LIB=/usr \
    RiVLib_DIR=/opt/rivlib/rivlib-2_6_0-x86_64-linux-gcc11

RUN mkdir -p /src/bin

# gedisimulator
RUN cd /src/gedisimulator && \
    for t in gediRat gediMetric mapLidar collocateWaves lasPoints \
             addNoiseHDF lgw2hdf fitTXpulse; do \
        make THIS="$t" && cp "$t" /src/bin/ || true; \
    done

# salcamake
RUN cd /src/salcamake && make THIS=hedgehog && cp hedgehog /src/bin/ || true
RUN cd /src/salcamake/cal_salca && make THIS=salcaCal && cp salcaCal /src/bin/ || true
RUN cd /src/salcamake/tiffSALCA && make THIS=tiffSALCA && cp tiffSALCA /src/bin/ || true

# voxelate
RUN cd /src/voxelate && \
    for t in voxelate voxelTLS voxDecimate dumpTiffAscii dumpTLScloud; do \
        make THIS="$t" && cp "$t" /src/bin/ || true; \
    done

# voxel_lidar
RUN cd /src/voxel_lidar && \
    for t in compareWaves optimiseTLS; do \
        make THIS="$t" && cp "$t" /src/bin/ || true; \
    done

# rxptobinary (CMake — clean build directory to avoid stale cache)
RUN rm -rf /src/rxptobinary/build && \
    mkdir -p /src/rxptobinary/build && \
    cmake -S /src/rxptobinary -B /src/rxptobinary/build \
          -DRiVLib_DIR="$RiVLib_DIR" && \
    cmake --build /src/rxptobinary/build && \
    cp /src/rxptobinary/build/readRXP /src/bin/ || true

# ── Runtime image ───────────────────────────────────────────────────────────
FROM fedora:43

RUN dnf install -y \
        gsl hdf5 gdal libgeotiff libtiff proj \
        python3 python3-numpy python3-h5py \
    && dnf clean all

# Riegl shared libraries
COPY --from=build /opt/rivlib/rivlib-2_6_0-x86_64-linux-gcc11/lib /opt/rivlib/lib
COPY --from=build /opt/rdblib/library /opt/rdblib/library
RUN echo "/opt/rivlib/lib" > /etc/ld.so.conf.d/rivlib.conf && \
    echo "/opt/rdblib/library" >> /etc/ld.so.conf.d/rivlib.conf && \
    ldconfig

# Binaries
COPY --from=build /src/bin/ /usr/local/bin/

# Python scripts
COPY voxelate/combinePAI.py /usr/local/bin/
COPY voxelate/changePAI.py  /usr/local/bin/

WORKDIR /data
ENTRYPOINT ["/bin/bash"]
