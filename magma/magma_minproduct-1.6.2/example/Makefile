# This is a simple standalone example. See README.txt
# Initially it is setup to use OpenBLAS.
# See magma_minproduct/make.inc for alternate BLAS and LAPACK libraries,
# or use pkg-config as described below.

# Paths where MAGMA_minproduct, CUDA, and OpenBLAS are installed
MAGMA_minproductDIR     ?= ..
CUDADIR      ?= /usr/local/cuda
OPENBLASDIR  ?= /usr/local/openblas

CC            = gcc
LD            = gcc
CFLAGS        = -Wall
LDFLAGS       = -Wall

MAGMA_minproduct_CFLAGS := -DADD_ -I$(MAGMA_minproductDIR)/include -I$(CUDADIR)/include

MAGMA_minproduct_LIBS   := -L$(MAGMA_minproductDIR)/lib -L$(CUDADIR)/lib64 -L$(OPENBLASDIR)/lib \
                -lmagma_minproduct -lcublas -lcudart -lopenblas

# Alternatively, using pkg-config (see README.txt):
# MAGMA_minproduct_CFLAGS := $(shell pkg-config --cflags magma_minproduct)
# MAGMA_minproduct_LIBS   := $(shell pkg-config --libs   magma_minproduct)

all: example

%.o: %.c
	$(CC) $(CFLAGS) $(MAGMA_minproduct_CFLAGS) -c -o $@ $<

example: example.o zfill.o
	$(LD) $(LDFLAGS) -o $@ $^ $(MAGMA_minproduct_LIBS)

clean:
	-rm -f example *.o
