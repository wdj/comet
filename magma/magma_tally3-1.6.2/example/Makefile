# This is a simple standalone example. See README.txt
# Initially it is setup to use OpenBLAS.
# See magma_tally3/make.inc for alternate BLAS and LAPACK libraries,
# or use pkg-config as described below.

# Paths where MAGMA_tally3, CUDA, and OpenBLAS are installed
MAGMA_tally3DIR     ?= ..
CUDADIR      ?= /usr/local/cuda
OPENBLASDIR  ?= /usr/local/openblas

CC            = gcc
LD            = gcc
CFLAGS        = -Wall
LDFLAGS       = -Wall

MAGMA_tally3_CFLAGS := -DADD_ -I$(MAGMA_tally3DIR)/include -I$(CUDADIR)/include

MAGMA_tally3_LIBS   := -L$(MAGMA_tally3DIR)/lib -L$(CUDADIR)/lib64 -L$(OPENBLASDIR)/lib \
                -lmagma_tally3 -lcublas -lcudart -lopenblas

# Alternatively, using pkg-config (see README.txt):
# MAGMA_tally3_CFLAGS := $(shell pkg-config --cflags magma_tally3)
# MAGMA_tally3_LIBS   := $(shell pkg-config --libs   magma_tally3)

all: example

%.o: %.c
	$(CC) $(CFLAGS) $(MAGMA_tally3_CFLAGS) -c -o $@ $<

example: example.o zfill.o
	$(LD) $(LDFLAGS) -o $@ $^ $(MAGMA_tally3_LIBS)

clean:
	-rm -f example *.o
