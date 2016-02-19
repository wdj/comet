#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA_tally3 (version 1.6.1) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      @date January 2015
#//////////////////////////////////////////////////////////////////////////////

MAGMA_tally3_DIR = ..
include ../Makefile.internal

ZHDR =  magma_tally3_z.h 	\
	magma_tally3_zc.h 	\
	magma_tally3blas_z.h	\
	magma_tally3blas_zc.h	\
	magma_tally3_zlapack.h

-include Makefile.src

ALLHDR= magma_tally3.h 	\
	magma_tally3blas.h	\
	magma_tally3_lapack.h	\
	$(ZHDR) $(CHDR) $(DHDR) $(SHDR)

all: $(ALLHDR)

clean:
	rm -f *~ *.mod

cleanall: clean

