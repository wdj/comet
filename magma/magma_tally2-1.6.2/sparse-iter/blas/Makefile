#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA_tally2 (version 1.6.2) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      @date May 2015
#//////////////////////////////////////////////////////////////////////////////

DIR = sparse-iter/blas
MAGMA_tally2_DIR = ../..
include $(MAGMA_tally2_DIR)/Makefile.internal

ZHDR := 

# alphabetic order by base name (ignoring precision)
ZSRC := \
	magma_tally2_z_blaswrapper.cpp	\
	zbajac_csr.cu			\
	zgecsrmv.cu				\
	zgeellmv.cu				\
	zgeelltmv.cu			\
	zgeellrtmv.cu			\
	zgesellcmv.cu			\
	zgesellcmmv.cu			\
	zjacobisetup.cu			\
	zlobpcg_shift.cu		\
	zlobpcg_residuals.cu	\
	zlobpcg_maxpy.cu		\
	zmdot.cu				\
	zmergebicgstab.cu		\
	zmergebicgstab2.cu		\
	zmergecg.cu				\
	zmgecsrmv.cu			\
	zmgeellmv.cu			\
	zmgeelltmv.cu			\
	zmgesellcmmv.cu			\
	zpipelinedgmres.cu		\


# Iterative factorizations
ZSRC += \
	ziterilu.cu	\
	ziteric.cu	\


# Auxiliary routines
ZSRC += \
	zcompact.cu				\
	magma_tally2_zmcsrcompressor_gpu.cu	\


# ------------------------------------------------------------
-include Makefile.local
-include Makefile.src

ALLHDR := $(ZHDR) $(CHDR) $(DHDR) $(SHDR)

ALLSRC := $(ZSRC) $(CSRC) $(DSRC) $(SSRC)

ALLOBJ := $(ALLSRC:.cu=.$(o_ext)) 
ALLOBJ := $(ALLOBJ:.cpp=.$(o_ext))

.PHONY: all clean cleanall

all: $(LIBMAGMA_tally2_SPARSE) 

$(LIBMAGMA_tally2_SPARSE): $(ALLHDR) $(ALLOBJ) stamp
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

# this forces LIBMAGMA_tally2 to be archived
stamp: $(ALLHDR) $(ALLOBJ)
	touch $@

clean:
	rm -f *.$(o_ext) *~ *.linkinfo

cleanall: clean
	rm -f $(LIBMAGMA_tally2_SPARSE)
