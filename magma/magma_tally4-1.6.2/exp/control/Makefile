#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA_tally4 (version 1.6.1) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      @date January 2015
#//////////////////////////////////////////////////////////////////////////////

MAGMA_tally4_DIR = ..
QUARKDIR = ../quark
include ../Makefile.internal

ifeq ($(PGIFORTRAN), 0)
ZSRC := $(ZSRC) 		\
	magma_tally4_zf77.cpp		\
	magma_tally4_param.f90		\
	magma_tally4_zfortran.f90	\
	magma_tally4.f90
else
ZSRC := $(ZSRC) 		\
	magma_tally4_zf77pgi.cpp
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

all: $(LIBMAGMA_tally4)

$(LIBMAGMA_tally4): $(ALLOBJ)
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
	rm -f $(LIBMAGMA_tally4)

# Extra dependencies
magma_tally4.o: magma_tally4_zfortran.o magma_tally4_cfortran.o magma_tally4_dfortran.o magma_tally4_sfortran.o 
