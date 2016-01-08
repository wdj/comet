#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA_tally4 (version 1.6.1) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      September 2013
#//////////////////////////////////////////////////////////////////////////////

MAGMA_tally4_DIR = ../..
include $(MAGMA_tally4_DIR)/Makefile.internal
include ../make.exp.inc

ZSRC := \

-include Makefile.local
-include Makefile.src

ALLSRC := $(ZSRC) $(CSRC) $(DSRC) $(SSRC) \
	ca_dbg_tools.cpp \
	Llog.cpp
	
ALLOBJ := $(ALLSRC:.cpp=.o)
ALLOBJ := $(ALLOBJ:.f90=.o)
ALLOBJ := $(ALLOBJ:.F90=.o)

.PHONY: all clean cleanall

all: $(LIBMAGMA_tally4_QUARK)

$(LIBMAGMA_tally4_QUARK): stamp
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

# this forces LIBMAGMA_tally4 to be archived
stamp: $(ALLOBJ)
	touch $@

clean:
	rm -f *.cu_o *.o *~ *.linkinfo stamp

cleanall: clean
	rm -f $(LIBMAGMA_tally4_QUARK)