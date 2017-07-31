#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA_tally2 (version 1.6.2) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      @date May 2015
#//////////////////////////////////////////////////////////////////////////////

DIR = sparse-iter/control
MAGMA_tally2_DIR = ../..
include $(MAGMA_tally2_DIR)/Makefile.internal

INC += -I../include

ZSRC := \
	error.cpp				\
	magma_tally2_zdomainoverlap.cpp    \
	magma_tally2_zutil_sparse.cpp	\
	magma_tally2_zfree.cpp		    \
	magma_tally2_zmatrixchar.cpp	\
	magma_tally2_zmconverter.cpp	\
	magma_tally2_zmgenerator.cpp	\
	magma_tally2_zmio.cpp			\
	magma_tally2_zsolverinfo.cpp	\
	magma_tally2_zcsrsplit.cpp		\
	magma_tally2_ziteriluutils.cpp	\
	magma_tally2_zmcsrpass.cpp		\
	magma_tally2_zmcsrpass_gpu.cpp	\
	magma_tally2_zmcsrcompressor.cpp	\
	magma_tally2_zmscale.cpp		\
	magma_tally2_zmdiff.cpp		\
	magma_tally2_zmlumerge.cpp		\
	magma_tally2_zmtranspose.cpp	\
	magma_tally2_zmtransfer.cpp	\
	magma_tally2_zmilustruct.cpp	\
	magma_tally2_zsort.cpp		\
	magma_tally2_zvinit.cpp		\
	magma_tally2_zvio.cpp			\
	magma_tally2_zvtranspose.cpp	\
	magma_tally2_zvpass.cpp		\
	magma_tally2_zvpass_gpu.cpp	\

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

all: $(LIBMAGMA_tally2_SPARSE)

$(LIBMAGMA_tally2_SPARSE): $(ALLOBJ) stamp
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

# this forces LIBMAGMA_tally2 to be archived
stamp: $(ALLHDR) $(ALLOBJ)
	touch $@

clean:
	rm -f *.$(o_ext) *~ *.linkinfo *.mod ../include/*.mod sizeptr

cleanall: clean
	rm -f $(LIBMAGMA_tally2_SPARSE)
