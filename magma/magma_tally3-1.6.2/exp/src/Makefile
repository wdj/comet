#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA_tally3 (version 1.6.1) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      @date January 2015
#//////////////////////////////////////////////////////////////////////////////

MAGMA_tally3_DIR = ..
QUARKDIR = ../../quark
LIBMAGMA_tally3BLAS = ../../lib/libmagma_tally3blas.a
include ../Makefile.internal

INC += -I$(MAGMA_tally3_DIR)/control

ZSRC  = zpotrf_mc.cpp	\
        zgetrf_mc.cpp   \
        zgeqrf_mc.cpp   \
        zgeqrf-v2.cpp   \
        zgeqrf-v3.cpp   \
        zlarfb_gpu.cpp   

-include Makefile.src

ALLSRC= $(ZSRC) $(CSRC) $(DSRC) $(SSRC)
ALLOBJ  = $(ALLSRC:.cpp=.o)

all: $(LIBMAGMA_tally3)

$(LIBMAGMA_tally3): $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

.cpp.o:
	$(CC) $(OPTS) $(INC) -c $< -o $@

clean:
	rm -f *.o *~ *.linkinfo

cleanall: clean
	rm -f *.a
	rm -f $(LIBMAGMA_tally3)
