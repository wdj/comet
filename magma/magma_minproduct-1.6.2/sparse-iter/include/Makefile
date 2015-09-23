#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA_minproduct (version 1.1) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      @date
#//////////////////////////////////////////////////////////////////////////////

DIR = sparse-iter/include
MAGMA_minproduct_DIR = ../..
include $(MAGMA_minproduct_DIR)/Makefile.internal

ZHDR := \
	magma_minproductsparse_z.h		\
	magma_minproductsparse_zc.h	\
	magma_minproductsparse_types.h	\


# ------------------------------------------------------------
-include Makefile.local
-include Makefile.src

ALLHDR := $(ZHDR) $(CHDR) $(DHDR) $(SHDR)

.PHONY: all clean cleanall

all: $(ALLHDR)

clean:
	rm -f *~ *.mod

cleanall: clean

