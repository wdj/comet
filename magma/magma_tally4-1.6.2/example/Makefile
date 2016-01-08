# This is a simple standalone example. See README.txt
# Initially it is setup to use OpenBLAS.
# See magma_tally4/make.inc for alternate BLAS and LAPACK libraries,
# or use pkg-config as described below.

# Paths where MAGMA_tally4, CUDA, and OpenBLAS are installed
MAGMA_tally4DIR     ?= ..
CUDADIR      ?= /usr/local/cuda
OPENBLASDIR  ?= /usr/local/openblas

CC            = gcc
LD            = gcc
CFLAGS        = -Wall
LDFLAGS       = -Wall

MAGMA_tally4_CFLAGS := -DADD_ -I$(MAGMA_tally4DIR)/include -I$(CUDADIR)/include

MAGMA_tally4_LIBS   := -L$(MAGMA_tally4DIR)/lib -L$(CUDADIR)/lib64 -L$(OPENBLASDIR)/lib \
                -lmagma_tally4 -lcublas -lcudart -lopenblas

# Alternatively, using pkg-config (see README.txt):
# MAGMA_tally4_CFLAGS := $(shell pkg-config --cflags magma_tally4)
# MAGMA_tally4_LIBS   := $(shell pkg-config --libs   magma_tally4)

all: example

%.o: %.c
	$(CC) $(CFLAGS) $(MAGMA_tally4_CFLAGS) -c -o $@ $<

example: example.o zfill.o
	$(LD) $(LDFLAGS) -o $@ $^ $(MAGMA_tally4_LIBS)

clean:
	-rm -f example *.o
