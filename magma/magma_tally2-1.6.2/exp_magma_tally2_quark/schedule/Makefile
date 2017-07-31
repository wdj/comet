#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA_tally2 (version 1.6.1) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      September 2013
#//////////////////////////////////////////////////////////////////////////////

MAGMA_tally2_DIR = ../..
include $(MAGMA_tally2_DIR)/Makefile.internal
include ../make.exp.inc

ZSRC := schedule.cpp
	
-include Makefile.local
-include Makefile.src

ALLSRC := $(ZSRC) $(CSRC) $(DSRC) $(SSRC)
		
ALLOBJ := $(ALLSRC:.cpp=.o)
ALLOBJ := $(ALLOBJ:.f90=.o)
ALLOBJ := $(ALLOBJ:.F90=.o)

.PHONY: all clean cleanall

all: $(LIBMAGMA_tally2_QUARK)

$(LIBMAGMA_tally2_QUARK): stamp
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

# this forces LIBMAGMA_tally2 to be archived
stamp: $(ALLOBJ)
	touch $@

clean:
	rm -f *.cu_o *.o *~ *.linkinfo stamp

cleanall: clean
	rm -f $(LIBMAGMA_tally2_QUARK)