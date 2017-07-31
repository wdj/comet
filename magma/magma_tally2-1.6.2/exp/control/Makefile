#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA_tally2 (version 1.6.1) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      @date January 2015
#//////////////////////////////////////////////////////////////////////////////

MAGMA_tally2_DIR = ..
QUARKDIR = ../quark
include ../Makefile.internal

ifeq ($(PGIFORTRAN), 0)
ZSRC := $(ZSRC) 		\
	magma_tally2_zf77.cpp		\
	magma_tally2_param.f90		\
	magma_tally2_zfortran.f90	\
	magma_tally2.f90
else
ZSRC := $(ZSRC) 		\
	magma_tally2_zf77pgi.cpp
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

all: $(LIBMAGMA_tally2)

$(LIBMAGMA_tally2): $(ALLOBJ)
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
	rm -f $(LIBMAGMA_tally2)

# Extra dependencies
magma_tally2.o: magma_tally2_zfortran.o magma_tally2_cfortran.o magma_tally2_dfortran.o magma_tally2_sfortran.o 
