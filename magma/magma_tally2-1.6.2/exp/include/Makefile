#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA_tally2 (version 1.6.1) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      @date January 2015
#//////////////////////////////////////////////////////////////////////////////

MAGMA_tally2_DIR = ..
include ../Makefile.internal

ZHDR =  magma_tally2_z.h 	\
	magma_tally2_zc.h 	\
	magma_tally2blas_z.h	\
	magma_tally2blas_zc.h	\
	magma_tally2_zlapack.h

-include Makefile.src

ALLHDR= magma_tally2.h 	\
	magma_tally2blas.h	\
	magma_tally2_lapack.h	\
	$(ZHDR) $(CHDR) $(DHDR) $(SHDR)

all: $(ALLHDR)

clean:
	rm -f *~ *.mod

cleanall: clean

