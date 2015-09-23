#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA_minproduct (version 1.6.2) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      @date May 2015
#//////////////////////////////////////////////////////////////////////////////

DIR = sparse-iter/control
MAGMA_minproduct_DIR = ../..
include $(MAGMA_minproduct_DIR)/Makefile.internal

INC += -I../include

ZSRC := \
	error.cpp				\
	magma_minproduct_zdomainoverlap.cpp    \
	magma_minproduct_zutil_sparse.cpp	\
	magma_minproduct_zfree.cpp		    \
	magma_minproduct_zmatrixchar.cpp	\
	magma_minproduct_zmconverter.cpp	\
	magma_minproduct_zmgenerator.cpp	\
	magma_minproduct_zmio.cpp			\
	magma_minproduct_zsolverinfo.cpp	\
	magma_minproduct_zcsrsplit.cpp		\
	magma_minproduct_ziteriluutils.cpp	\
	magma_minproduct_zmcsrpass.cpp		\
	magma_minproduct_zmcsrpass_gpu.cpp	\
	magma_minproduct_zmcsrcompressor.cpp	\
	magma_minproduct_zmscale.cpp		\
	magma_minproduct_zmdiff.cpp		\
	magma_minproduct_zmlumerge.cpp		\
	magma_minproduct_zmtranspose.cpp	\
	magma_minproduct_zmtransfer.cpp	\
	magma_minproduct_zmilustruct.cpp	\
	magma_minproduct_zsort.cpp		\
	magma_minproduct_zvinit.cpp		\
	magma_minproduct_zvio.cpp			\
	magma_minproduct_zvtranspose.cpp	\
	magma_minproduct_zvpass.cpp		\
	magma_minproduct_zvpass_gpu.cpp	\

SRC := \
	mmio.cpp				\


# ------------------------------------------------------------
-include Makefile.local
-include Makefile.src

ALLSRC := $(ZSRC) $(CSRC) $(DSRC) $(SSRC) $(SRC)

ALLOBJ := $(ALLSRC:.cpp=.$(o_ext))
ALLOBJ := $(ALLOBJ:.f90=.$(o_ext))
ALLOBJ := $(ALLOBJ:.F90=.$(o_ext))

.PHONY: all clean cleanall

all: $(LIBMAGMA_minproduct_SPARSE)

$(LIBMAGMA_minproduct_SPARSE): $(ALLOBJ) stamp
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

# this forces LIBMAGMA_minproduct to be archived
stamp: $(ALLHDR) $(ALLOBJ)
	touch $@

clean:
	rm -f *.$(o_ext) *~ *.linkinfo *.mod ../include/*.mod sizeptr

cleanall: clean
	rm -f $(LIBMAGMA_minproduct_SPARSE)
