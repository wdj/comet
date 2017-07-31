# This is a simple standalone example. See README.txt
# Initially it is setup to use OpenBLAS.
# See magma_tally2/make.inc for alternate BLAS and LAPACK libraries,
# or use pkg-config as described below.

# Paths where MAGMA_tally2, CUDA, and OpenBLAS are installed
MAGMA_tally2DIR     ?= ..
CUDADIR      ?= /usr/local/cuda
OPENBLASDIR  ?= /usr/local/openblas

CC            = gcc
LD            = gcc
CFLAGS        = -Wall
LDFLAGS       = -Wall

MAGMA_tally2_CFLAGS := -DADD_ -I$(MAGMA_tally2DIR)/include -I$(CUDADIR)/include

MAGMA_tally2_LIBS   := -L$(MAGMA_tally2DIR)/lib -L$(CUDADIR)/lib64 -L$(OPENBLASDIR)/lib \
                -lmagma_tally2 -lcublas -lcudart -lopenblas

# Alternatively, using pkg-config (see README.txt):
# MAGMA_tally2_CFLAGS := $(shell pkg-config --cflags magma_tally2)
# MAGMA_tally2_LIBS   := $(shell pkg-config --libs   magma_tally2)

all: example

%.o: %.c
	$(CC) $(CFLAGS) $(MAGMA_tally2_CFLAGS) -c -o $@ $<

example: example.o zfill.o
	$(LD) $(LDFLAGS) -o $@ $^ $(MAGMA_tally2_LIBS)

clean:
	-rm -f example *.o
