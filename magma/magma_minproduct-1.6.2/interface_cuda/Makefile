#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA_minproduct (version 1.6.1) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      @date January 2015
#//////////////////////////////////////////////////////////////////////////////

DIR = interface_cuda
MAGMA_minproduct_DIR = ..
include $(MAGMA_minproduct_DIR)/Makefile.internal

# alphabetic order by base name (ignoring precision)
ZSRC := \
	blas_z.cpp	\
	zset_get.cpp	\

SRC := \
	alloc.cpp	\
	error.cpp	\
	interface.cpp	\
	set_get.cpp	\
	error.cpp	\

-include Makefile.local
-include Makefile.src

ALLSRC := $(ZSRC) $(CSRC) $(DSRC) $(SSRC) $(SRC)

ALLOBJ := $(ALLSRC:.cpp=.o)

.PHONY: all clean cleanall

.DEFAULT_GOAL := all
all: $(LIBMAGMA_minproduct)

$(LIBMAGMA_minproduct): $(ALLOBJ) stamp
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

# this forces LIBMAGMA_minproduct to be archived
stamp: $(ALLOBJ)
	touch $@

clean:
	rm -f *.o *~ *.linkinfo stamp

cleanall: clean
	rm -f $(LIBMAGMA_minproduct)
