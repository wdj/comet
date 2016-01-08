#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA_tally4 (version 1.1) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      @date
#//////////////////////////////////////////////////////////////////////////////

DIR = sparse-iter/include
MAGMA_tally4_DIR = ../..
include $(MAGMA_tally4_DIR)/Makefile.internal

ZHDR := \
	magma_tally4sparse_z.h		\
	magma_tally4sparse_zc.h	\
	magma_tally4sparse_types.h	\


# ------------------------------------------------------------
-include Makefile.local
-include Makefile.src

ALLHDR := $(ZHDR) $(CHDR) $(DHDR) $(SHDR)

.PHONY: all clean cleanall

all: $(ALLHDR)

clean:
	rm -f *~ *.mod

cleanall: clean

