#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA_minproduct (version 1.6.1) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      @date January 2015
#//////////////////////////////////////////////////////////////////////////////

MAGMA_minproduct_DIR = ..
include ../Makefile.internal

ZHDR =  magma_minproduct_z.h 	\
	magma_minproduct_zc.h 	\
	magma_minproductblas_z.h	\
	magma_minproductblas_zc.h	\
	magma_minproduct_zlapack.h

-include Makefile.src

ALLHDR= magma_minproduct.h 	\
	magma_minproductblas.h	\
	magma_minproduct_lapack.h	\
	$(ZHDR) $(CHDR) $(DHDR) $(SHDR)

all: $(ALLHDR)

clean:
	rm -f *~ *.mod

cleanall: clean

