FROM fedora:43

# ── System packages ──────────────────────────────────────────────────────────
RUN dnf install -y \
        gcc gcc-c++ make cmake \
        gsl gsl-devel \
        hdf5 hdf5-devel \
        gdal gdal-devel \
        libgeotiff libgeotiff-devel \
        libtiff libtiff-devel \
        proj proj-devel \
        python3 python3-pip python3-numpy python3-h5py \
        git wget unzip tcsh \
    && dnf clean all

# ── Clone external dependencies ──────────────────────────────────────────────
# Hancock C-tools library
RUN git clone https://bitbucket.org/StevenHancock/tools.git /opt/tools

# cmpfit / Minpack (Levenberg-Marquardt)
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

# Ensure the rivlib shared libs are discoverable
RUN echo "/opt/rivlib/rivlib-2_6_0-x86_64-linux-gcc11/lib" > /etc/ld.so.conf.d/rivlib.conf && \
    echo "/opt/rdblib/lib"                                 >> /etc/ld.so.conf.d/rivlib.conf && \
    ldconfig

# ── Environment variables ────────────────────────────────────────────────────
ENV ARCH=x86_64 \
    HANCOCKTOOLS_ROOT=/opt/tools \
    CMPFIT_ROOT=/opt/cmpfit \
    LIBCLIDAR_ROOT=/workspaces/libClidar/libclidar \
    VOXELATE_ROOT=/workspaces/libClidar/voxelate \
    GEDIRAT_ROOT=/workspaces/libClidar/gedisimulator \
    GSL_ROOT=/usr/lib64 \
    HDF5_LIB=/usr \
    RiVLib_DIR=/opt/rivlib/rivlib-2_6_0-x86_64-linux-gcc11

ENV PATH="${PATH}:/usr/local/bin:/workspaces/libClidar/bin"

# ── Output bin directory ─────────────────────────────────────────────────────
RUN mkdir -p /workspaces/libClidar/bin

WORKDIR /workspaces/libClidar
