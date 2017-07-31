#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA_tally2 (version 1.6.1) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      @date January 2015
#//////////////////////////////////////////////////////////////////////////////

MAGMA_tally2_DIR = ..
QUARKDIR = ../quark
LIBMAGMA_tally2BLAS = ../../lib/libmagma_tally2blas.a
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

$(EXE): $(LIBMAGMA_tally2) $(LIBMAGMA_tally2BLAS)

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
# MAGMA_tally2 Fortran module. Once finalized, the generated magma_tally2.mod
# has to be moved in the include directory
#
testing_%_f : testing_%_f.o $(FOBJ)
	$(FORT) $(OPTS) $(LDOPTS) $< $(FOBJ) -o $@ -L../lib -L../../lib \
          -lcuda -lmagma_tally2 -lmagma_tally2blas -lmagma_tally2 $(LIBDIR) $(LIB) 

testing_% : testing_%.o
	$(CC) $(OPTS) $(LDOPTS) $(NOMAIN) $< $(OBJ) -o $@ -L../lib -L../../lib \
          -lcuda -lmagma_tally2 -lmagma_tally2blas -lmagma_tally2 $(LIBDIR) $(LIB) 

#
# These two files are provided by cuda
#
fortran.o : $(CUDADIR)/src/fortran.c
	$(CC) $(OPTS) $(INC) -c $< -o $@

fortran_thunking.o : $(CUDADIR)/src/fortran_thunking.c
	$(CC) $(OPTS) $(INC) -c $< -o $@
