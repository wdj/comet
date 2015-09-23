#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA_minproduct (version 1.6.1) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      @date January 2015
#//////////////////////////////////////////////////////////////////////////////

MAGMA_minproduct_DIR = ..
QUARKDIR = ../quark
LIBMAGMA_minproductBLAS = ../../lib/libmagma_minproductblas.a
include ../Makefile.internal

ZSRC  = testing_zpotrf_mc.cpp \
        testing_zgetrf_mc.cpp \
        testing_zgeqrf_mc.cpp \
        testing_zgeqrf-v2.cpp 
        

#testing_zswap.cpp

ZSRC += $(ZSRCF)
-include Makefile.src

ALLSRC  = $(ZSRC) $(CSRC) $(DSRC) $(SSRC)
ALLOBJF  = $(ALLSRC:.f90=.o)
ALLOBJF := $(ALLOBJF:.cuf=.o)
ALLOBJ   = $(ALLOBJF:.cpp=.o)
EXE      = $(ALLOBJ:.o=)

all: $(EXE)

$(EXE): $(LIBMAGMA_minproduct) $(LIBMAGMA_minproductBLAS)

clean:
	rm -f *.o *~

cleanall: clean
	rm -f $(EXE)

ifeq ($(FORT), pgfortran)
testing_%.o : testing_%.cuf
	$(FORT) $(FOPTS) $(INC) -c $< -o $@ 
else
testing_%.o : testing_%.f90
	$(FORT) $(FOPTS) $(INC) -c $< -o $@ 
endif

%.o: %.cpp
	$(CC) $(OPTS) $(INC) -c $< -o $@

#
# MAGMA_minproduct Fortran module. Once finalized, the generated magma_minproduct.mod
# has to be moved in the include directory
#
testing_%_f : testing_%_f.o $(FOBJ)
	$(FORT) $(OPTS) $(LDOPTS) $< $(FOBJ) -o $@ -L../lib -L../../lib \
          -lcuda -lmagma_minproduct -lmagma_minproductblas -lmagma_minproduct $(LIBDIR) $(LIB) 

testing_% : testing_%.o
	$(CC) $(OPTS) $(LDOPTS) $(NOMAIN) $< $(OBJ) -o $@ -L../lib -L../../lib \
          -lcuda -lmagma_minproduct -lmagma_minproductblas -lmagma_minproduct $(LIBDIR) $(LIB) 

#
# These two files are provided by cuda
#
fortran.o : $(CUDADIR)/src/fortran.c
	$(CC) $(OPTS) $(INC) -c $< -o $@

fortran_thunking.o : $(CUDADIR)/src/fortran_thunking.c
	$(CC) $(OPTS) $(INC) -c $< -o $@
