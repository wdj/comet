#==============================================================================
# Makefile for genomics metric calculation code.

SHELL   = /bin/sh
#CC      = gcc
CC      = cc

#CFLAGS  += -pedantic -Wall -O2
CFLAGS  += -Wall -Wno-unused-function -O2 -finline-functions -finline-limit=1000
CFLAGS  += -fargument-noalias-anything -O2
CFLAGS  += -I../magma/magma_minproduct-1.6.2/include $(CRAY_CUDATOOLKIT_INCLUDE_OPTS) -DADD_
CFLAGS  += -I../src -I.
CFLAGS += -DFP_PRECISION_DOUBLE
#CFLAGS += -DFP_PRECISION_SINGLE

LFLAGS = $(LDFLAGS) -L../magma/magma_minproduct-1.6.2/lib -lmagma_minproduct $(CRAY_CUDATOOLKIT_POST_LINK_OPTS) -lcublas

TARGET  = genomics_metric
SOURCES = $(shell echo ../src/*.c ../driver/*.c)
HEADERS = $(shell echo ../src/*.h ../driver/*.h)
OBJECTS = $(SOURCES:.c=.o)

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

%.o: %.c $(HEADERS)
# %.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	-rm -f $(OBJECTS)
 
distclean: clean
	-rm -f $(TARGET)

.PHONY: all clean distclean
#==============================================================================
