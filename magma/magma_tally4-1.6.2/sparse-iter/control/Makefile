#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA_tally4 (version 1.6.2) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      @date May 2015
#//////////////////////////////////////////////////////////////////////////////

DIR = sparse-iter/control
MAGMA_tally4_DIR = ../..
include $(MAGMA_tally4_DIR)/Makefile.internal

INC += -I../include

ZSRC := \
	error.cpp				\
	magma_tally4_zdomainoverlap.cpp    \
	magma_tally4_zutil_sparse.cpp	\
	magma_tally4_zfree.cpp		    \
	magma_tally4_zmatrixchar.cpp	\
	magma_tally4_zmconverter.cpp	\
	magma_tally4_zmgenerator.cpp	\
	magma_tally4_zmio.cpp			\
	magma_tally4_zsolverinfo.cpp	\
	magma_tally4_zcsrsplit.cpp		\
	magma_tally4_ziteriluutils.cpp	\
	magma_tally4_zmcsrpass.cpp		\
	magma_tally4_zmcsrpass_gpu.cpp	\
	magma_tally4_zmcsrcompressor.cpp	\
	magma_tally4_zmscale.cpp		\
	magma_tally4_zmdiff.cpp		\
	magma_tally4_zmlumerge.cpp		\
	magma_tally4_zmtranspose.cpp	\
	magma_tally4_zmtransfer.cpp	\
	magma_tally4_zmilustruct.cpp	\
	magma_tally4_zsort.cpp		\
	magma_tally4_zvinit.cpp		\
	magma_tally4_zvio.cpp			\
	magma_tally4_zvtranspose.cpp	\
	magma_tally4_zvpass.cpp		\
	magma_tally4_zvpass_gpu.cpp	\

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

all: $(LIBMAGMA_tally4_SPARSE)

$(LIBMAGMA_tally4_SPARSE): $(ALLOBJ) stamp
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

# this forces LIBMAGMA_tally4 to be archived
stamp: $(ALLHDR) $(ALLOBJ)
	touch $@

clean:
	rm -f *.$(o_ext) *~ *.linkinfo *.mod ../include/*.mod sizeptr

cleanall: clean
	rm -f $(LIBMAGMA_tally4_SPARSE)
