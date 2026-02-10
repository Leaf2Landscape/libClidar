# Makefile for iter_thresh
#
HOME = /Users/dill
LIBS = -lm -lgsl -lgslcblas
VOXDIR=${VOXELATE_ROOT}
INCLS = -I/usr/local/include -I${HANCOCKTOOLS_ROOT} -I${VOXDIR} -I/usr/include/geotiff -I${CMPFIT_ROOT} -I${LIBCLIDAR_ROOT}
CFLAGS += -Wall
#CFLAGS += -g
CFLAGS += -O3
LIBFILES = $(LIBCLIDAR_ROOT)/gaussFit.o $(LIBCLIDAR_ROOT)/libLasProcess.o $(LIBCLIDAR_ROOT)/libLasRead.o $(LIBCLIDAR_ROOT)/libLidVoxel.o $(LIBCLIDAR_ROOT)/libTLSread.o
LOCLIB = gaussFit.o libLasProcess.o libLasRead.o libLidVoxel.o libTLSread.o
VOXLIST=$(VOXDIR)/voxFunc.o
VOX=$(subst $(VOXDIR)/,,$(VOXLIST))
LOCLIST=libOptWaveform.o

MINBITS=mpfit.o
MIN=$(CMPFIT_ROOT)/$(MINBITS)


CC = gcc
#CC= /opt/SUNWspro/bin/cc

#CFLAGS += -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64

THIS=compareWaves

$(THIS):	$(THIS).o $(VOXLIST) $(MIN) $(LOCLIST) $(LIBFILES)
		$(CC) $(CFLAGS) $@.o $(VOX) $(MINBITS) $(LOCLIST) $(LOCLIB) -o $@ $(LIBS) $(CFLAGS) $(INCLS)

.c.o:		$<
		$(CC) $(CFLAGS) -I. $(INCLS) -D$(ARCH)  -c $<

clean:
		\rm -f *% *~ *.o $(VOXDIR)/*.o

install:
		touch $(HOME)/bin/$(ARCH)/$(THIS)
		mv $(HOME)/bin/$(ARCH)/$(THIS) $(HOME)/bin/$(ARCH)/$(THIS).old
		cp $(THIS) $(HOME)/bin/$(ARCH)/$(THIS)
			
