# Makefile for iter_thresh
#
LIBS = -lm
INCLS = -I/usr/local/include -I${HANCOCKTOOLS_ROOT} -I${CMPFIT_ROOT}
CFLAGS += -Wall
CFLAGS += -g
#CFLAGS += -O3

MINBITS=mpfit.o
MIN=${CMPFIT_ROOT}/$(MINBITS)


CC = gcc
#CC= /opt/SUNWspro/bin/cc

#CFLAGS += -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64

THIS=salcaCal

$(THIS):	$(THIS).o $(MIN)
		$(CC) $(CFLAGS) $@.o $(MINBITS) -o $@ $(LIBS) $(CFLAGS) $(INCLS)

.c.o:		$<
		$(CC) $(CFLAGS) -I. $(INCLS) -D$(ARCH)  -c $<

clean:
		\rm -f *% *~ *.o

install:
		touch $(HOME)/bin/$(ARCH)/$(THIS)
		mv $(HOME)/bin/$(ARCH)/$(THIS) $(HOME)/bin/$(ARCH)/$(THIS).old
		cp $(THIS) $(HOME)/bin/$(ARCH)/$(THIS)
			
