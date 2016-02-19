#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA_tally3 (version 1.1) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      @date
#//////////////////////////////////////////////////////////////////////////////

DIR = sparse-iter/include
MAGMA_tally3_DIR = ../..
include $(MAGMA_tally3_DIR)/Makefile.internal

ZHDR := \
	magma_tally3sparse_z.h		\
	magma_tally3sparse_zc.h	\
	magma_tally3sparse_types.h	\


# ------------------------------------------------------------
-include Makefile.local
-include Makefile.src

ALLHDR := $(ZHDR) $(CHDR) $(DHDR) $(SHDR)

.PHONY: all clean cleanall

all: $(ALLHDR)

clean:
	rm -f *~ *.mod

cleanall: clean

