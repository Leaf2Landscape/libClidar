# Makefile for iter_thresh
#
LIBS = -lm -ltiff -lgeotiff -lgsl -lgslcblas -lhdf5 -L${HDF5_LIB}/lib #-lhdfeos -lGctp -lmfhdf -ldf -ljpeg -lz -lm -lsz 
MIN=mpfit.o
INCLS = -I/usr/local/include -I$(HANCOCKTOOLS_ROOT) -I${GSL_ROOT} -I${GSL_ROOT}/fft -I/usr/include/libgeotiff -I/opt/local/include -I$(LIBCLIDAR_ROOT) -I$(CMPFIT_ROOT) -I${HDF5_LIB}/include -I/usr/include/geotiff
CFLAGS += -Wall
#CFLAGS += -g
CFLAGS += -O3
LIBFILES = $(LIBCLIDAR_ROOT)/tiffWrite.o $(LIBCLIDAR_ROOT)/gaussFit.o $(LIBCLIDAR_ROOT)/libLasProcess.o $(LIBCLIDAR_ROOT)/libLasRead.o $(LIBCLIDAR_ROOT)/libLidVoxel.o $(LIBCLIDAR_ROOT)/libTLSread.o $(LIBCLIDAR_ROOT)/tiffRead.o $(LIBCLIDAR_ROOT)/libDEMhandle.o $(LIBCLIDAR_ROOT)/libLidarHDF.o
LOCLIB = tiffWrite.o gaussFit.o libLasProcess.o libLasRead.o libLidVoxel.o libTLSread.o tiffRead.o libDEMhandle.o libLidarHDF.o
GSLFit=linear.o 
MIN=mpfit.o

CC = gcc
#CC= /opt/SUNWspro/bin/cc

#CFLAGS += -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64

THIS=voxelate

$(THIS):	$(THIS).o ${CMPFIT_ROOT}/$(MIN) $(LIBFILES)
		$(CC) $(CFLAGS) $@.o $(LOCAL_OTHERS) $(MIN) $(LOCLIB) -o $@ $(LIBS) $(CFLAGS) $(INCLS)

.c.o:		$<
		$(CC) $(CFLAGS) -I. $(INCLS) -D$(ARCH)  -c $<

clean:
		\rm -f *% *~ *.o #$(TOOLDIR)/*.o

install:
		touch $(HOME)/bin/$(ARCH)/$(THIS)
		mv $(HOME)/bin/$(ARCH)/$(THIS) $(HOME)/bin/$(ARCH)/$(THIS).old
		cp $(THIS) $(HOME)/bin/$(ARCH)/$(THIS)

