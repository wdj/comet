#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA_minproduct (version 1.6.1) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      @date January 2015
#//////////////////////////////////////////////////////////////////////////////

MAGMA_minproduct_DIR = ..
QUARKDIR = ../../quark
LIBMAGMA_minproductBLAS = ../../lib/libmagma_minproductblas.a
include ../Makefile.internal

INC += -I$(MAGMA_minproduct_DIR)/control

ZSRC  = zpotrf_mc.cpp	\
        zgetrf_mc.cpp   \
        zgeqrf_mc.cpp   \
        zgeqrf-v2.cpp   \
        zgeqrf-v3.cpp   \
        zlarfb_gpu.cpp   

-include Makefile.src

ALLSRC= $(ZSRC) $(CSRC) $(DSRC) $(SSRC)
ALLOBJ  = $(ALLSRC:.cpp=.o)

all: $(LIBMAGMA_minproduct)

$(LIBMAGMA_minproduct): $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

.cpp.o:
	$(CC) $(OPTS) $(INC) -c $< -o $@

clean:
	rm -f *.o *~ *.linkinfo

cleanall: clean
	rm -f *.a
	rm -f $(LIBMAGMA_minproduct)
