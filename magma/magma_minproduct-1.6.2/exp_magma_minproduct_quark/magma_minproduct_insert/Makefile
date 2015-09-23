#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA_minproduct (version 1.6.1) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      September 2013
#//////////////////////////////////////////////////////////////////////////////

MAGMA_minproduct_DIR = ../..
include $(MAGMA_minproduct_DIR)/Makefile.internal
include ../make.exp.inc

ZSRC := magma_minproduct_insert_core_d.cpp \
	magma_minproduct_insert_d.cpp \
	magma_minproduct_insert_dev_d.cpp \

-include Makefile.local
-include Makefile.src

ALLSRC := $(ZSRC) $(CSRC) $(DSRC) $(SSRC)

ALLOBJ := $(ALLSRC:.cpp=.o)
ALLOBJ := $(ALLOBJ:.f90=.o)
ALLOBJ := $(ALLOBJ:.F90=.o)

.PHONY: all clean cleanall

all: $(LIBMAGMA_minproduct_QUARK)

$(LIBMAGMA_minproduct_QUARK): stamp
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

# this forces LIBMAGMA_minproduct to be archived
stamp: $(ALLOBJ)
	touch $@

clean:
	rm -f *.cu_o *.o *~ *.linkinfo stamp

cleanall: clean
	rm -f $(LIBMAGMA_minproduct_QUARK)
