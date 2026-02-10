# Makefile for iter_thresh
#
LIBS = -lm # -ltiff #-lhdfeos -lGctp -lmfhdf -ldf -ljpeg -lz -lm -lsz 
#NR_DIR=${HOME}/src/recipes
INCLS = -I/usr/local/include -I${HANCOCKTOOLS_ROOT} # -I${NR_DIR}
CFLAGS += -Wall
#CFLAGS += -g
CFLAGS += -O3
LOCAL_OTHERS =  #nrutil.o four1.o mrqminMine.o mrqcofMine.o gaussjMine.c covsrt.o


CC = gcc
#CC= /opt/SUNWspro/bin/cc

#CFLAGS += -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64

THIS=hedgehog

$(THIS):	$(THIS).o $(LOCAL_OTHERS) # $(NR) $(TOOLDIR)/$(TOOLS).o
		$(CC) $(CFLAGS) $@.o $(LOCAL_OTHERS) -o $@ $(LIBS) $(CFLAGS) $(INCLS) # $(NRhere)

.c.o:		$<
		$(CC) $(CFLAGS) -I. $(INCLS) -D$(ARCH)  -c $<

clean:
		\rm -f *% *~ *.o #$(TOOLDIR)/*.o

install:
		touch $(HOME)/bin/$(ARCH)/$(THIS)
		mv $(HOME)/bin/$(ARCH)/$(THIS) $(HOME)/bin/$(ARCH)/$(THIS).old
		cp $(THIS) $(HOME)/bin/$(ARCH)/$(THIS)


			
