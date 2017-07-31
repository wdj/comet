#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA_tally2 (version 1.1) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      @date
#//////////////////////////////////////////////////////////////////////////////

DIR = sparse-iter/include
MAGMA_tally2_DIR = ../..
include $(MAGMA_tally2_DIR)/Makefile.internal

ZHDR := \
	magma_tally2sparse_z.h		\
	magma_tally2sparse_zc.h	\
	magma_tally2sparse_types.h	\


# ------------------------------------------------------------
-include Makefile.local
-include Makefile.src

ALLHDR := $(ZHDR) $(CHDR) $(DHDR) $(SHDR)

.PHONY: all clean cleanall

all: $(ALLHDR)

clean:
	rm -f *~ *.mod

cleanall: clean

