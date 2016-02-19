#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA_tally3 (version 1.6.2) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      @date May 2015
#//////////////////////////////////////////////////////////////////////////////

DIR = sparse-iter/control
MAGMA_tally3_DIR = ../..
include $(MAGMA_tally3_DIR)/Makefile.internal

INC += -I../include

ZSRC := \
	error.cpp				\
	magma_tally3_zdomainoverlap.cpp    \
	magma_tally3_zutil_sparse.cpp	\
	magma_tally3_zfree.cpp		    \
	magma_tally3_zmatrixchar.cpp	\
	magma_tally3_zmconverter.cpp	\
	magma_tally3_zmgenerator.cpp	\
	magma_tally3_zmio.cpp			\
	magma_tally3_zsolverinfo.cpp	\
	magma_tally3_zcsrsplit.cpp		\
	magma_tally3_ziteriluutils.cpp	\
	magma_tally3_zmcsrpass.cpp		\
	magma_tally3_zmcsrpass_gpu.cpp	\
	magma_tally3_zmcsrcompressor.cpp	\
	magma_tally3_zmscale.cpp		\
	magma_tally3_zmdiff.cpp		\
	magma_tally3_zmlumerge.cpp		\
	magma_tally3_zmtranspose.cpp	\
	magma_tally3_zmtransfer.cpp	\
	magma_tally3_zmilustruct.cpp	\
	magma_tally3_zsort.cpp		\
	magma_tally3_zvinit.cpp		\
	magma_tally3_zvio.cpp			\
	magma_tally3_zvtranspose.cpp	\
	magma_tally3_zvpass.cpp		\
	magma_tally3_zvpass_gpu.cpp	\

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

all: $(LIBMAGMA_tally3_SPARSE)

$(LIBMAGMA_tally3_SPARSE): $(ALLOBJ) stamp
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

# this forces LIBMAGMA_tally3 to be archived
stamp: $(ALLHDR) $(ALLOBJ)
	touch $@

clean:
	rm -f *.$(o_ext) *~ *.linkinfo *.mod ../include/*.mod sizeptr

cleanall: clean
	rm -f $(LIBMAGMA_tally3_SPARSE)
