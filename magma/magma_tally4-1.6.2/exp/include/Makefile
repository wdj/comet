#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA_tally4 (version 1.6.1) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      @date January 2015
#//////////////////////////////////////////////////////////////////////////////

MAGMA_tally4_DIR = ..
include ../Makefile.internal

ZHDR =  magma_tally4_z.h 	\
	magma_tally4_zc.h 	\
	magma_tally4blas_z.h	\
	magma_tally4blas_zc.h	\
	magma_tally4_zlapack.h

-include Makefile.src

ALLHDR= magma_tally4.h 	\
	magma_tally4blas.h	\
	magma_tally4_lapack.h	\
	$(ZHDR) $(CHDR) $(DHDR) $(SHDR)

all: $(ALLHDR)

clean:
	rm -f *~ *.mod

cleanall: clean

