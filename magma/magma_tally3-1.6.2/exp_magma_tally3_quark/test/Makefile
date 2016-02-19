#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA_tally3 (version 1.6.1) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      September 2013
#//////////////////////////////////////////////////////////////////////////////

MAGMA_tally3_DIR = ../..
include $(MAGMA_tally3_DIR)/Makefile.internal
include ../make.exp.inc

ZSRC := \
     testing_dgetrf_gpu_amc.cpp \
     testing_dgetrf_gpu_amc_v2.cpp \
	 testing_dgetrf_mgpu_amc.cpp \
	 testing_dgetrf_mgpu_amc_v2.cpp \
	 testing_dgetrf_mgpu_amc_v3.cpp \
	 
-include Makefile.local
-include Makefile.src

ALLSRC := $(ZSRC) $(CSRC) $(DSRC) $(SSRC)
	
ALLOBJ := $(ALLSRC:.cpp=.o)
ALLOBJ := $(ALLOBJ:.f90=.o)
ALLOBJ := $(ALLOBJ:.F90=.o)
EXE    := $(ALLOBJ:.o=)

.PHONY: all clean cleanall

all: $(EXE)

clean:
	rm -f *.cu_o *.o *~ *.linkinfo stamp $(EXE)

cleanall: clean
	rm -f $(LIBMAGMA_tally3)

testing_% : testing_%.o
	$(CC) $(LDOPTS) $(NOMAIN) $< -o $@ \
	$(LIBTEST) $(LIBLAPACKTEST) -L../lib -lmagma_tally3 \
	$(LIBDIR) \
	$(LIB)
