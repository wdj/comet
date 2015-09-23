#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA_minproduct (version 1.6.1) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      @date January 2015
#//////////////////////////////////////////////////////////////////////////////

MAGMA_minproduct_DIR = ..
QUARKDIR = ../quark
include ../Makefile.internal

ifeq ($(PGIFORTRAN), 0)
ZSRC := $(ZSRC) 		\
	magma_minproduct_zf77.cpp		\
	magma_minproduct_param.f90		\
	magma_minproduct_zfortran.f90	\
	magma_minproduct.f90
else
ZSRC := $(ZSRC) 		\
	magma_minproduct_zf77pgi.cpp
endif

-include Makefile.src

ALLSRC= $(ZSRC) $(CSRC) $(DSRC) $(SSRC)	\
        auto_tune.cpp    \
        init.cpp        \
        auxiliary.cpp 	\
        xerbla.cpp 	\
	get_nb.cpp

ALLOBJ  = $(ALLSRC:.cpp=.o)
ALLOBJ := $(ALLOBJ:.f90=.o)

all: $(LIBMAGMA_minproduct)

$(LIBMAGMA_minproduct): $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@
	cp *.mod ../include

.cpp.o:
	$(CC) $(OPTS) $(INC) -c $< -o $@

%.o : %.f90 
	$(FORT) $(FOPTS) $(INC) -c $< -o $@

clean:
	rm -f *.cu_o *.o *~ *.linkinfo *.mod

cleanall: clean
	rm -f *.a
	rm -f $(LIBMAGMA_minproduct)

# Extra dependencies
magma_minproduct.o: magma_minproduct_zfortran.o magma_minproduct_cfortran.o magma_minproduct_dfortran.o magma_minproduct_sfortran.o 
