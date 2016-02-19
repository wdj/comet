#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA_tally3 (version 1.6.2) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      @date May 2015
#//////////////////////////////////////////////////////////////////////////////

DIR = sparse-iter/testing
MAGMA_tally3_DIR = ../..
include $(MAGMA_tally3_DIR)/Makefile.internal

LIB += -lcusparse 

# ----------
# utility functions
ZSRC += \
	testing_zblas.cpp			\
	testing_zmatrix.cpp			\
	testing_zio.cpp				\
	testing_zmcompressor.cpp	\
	testing_zmconverter.cpp		\
	testing_zsort.cpp	\



# ----------
# low level LA operations
ZSRC += \
	testing_zdot.cpp			\
	testing_zspmv.cpp			\
	testing_zspmm.cpp			\
	testing_zmadd.cpp			\

# ----------
# iterative solvers
ZSRC += \
	testing_zsolver.cpp			\


# ------------------------------------------------------------
-include Makefile.local
-include Makefile.src

# path for testings.h, flops.h
INC += -I$(MAGMA_tally3_DIR)/testing

ALLSRC := $(ZSRC) $(CSRC) $(DSRC) $(SSRC)
ALLOBJ := $(ALLSRC:.cpp=.$(o_ext))
EXE    := $(ALLOBJ:.$(o_ext)=)

LIBTEST := $(MAGMA_tally3_DIR)/testing/libtest.a

.PHONY: all lib clean cleanall

# --------------------
all: $(EXE)

echo:
	@echo EXE $(EXE)
	@echo ALLOBJ $(ALLOBJ)
	@echo ALLSRC $(ALLSRC)
	@echo ZSRC $(ZSRC)


# build libraries if they have not yet been built
$(LIBTEST):
	@echo ======================================== libtest.a
	cd ../../testing && $(MAKE) libtest.a

$(LIBMAGMA_tally3):
	@echo ======================================== libmagma_tally3.a
	cd ../.. && $(MAKE) lib

$(LIBMAGMA_tally3_SPARSE):
	@echo ======================================== libmagma_tally3_sparse.a
	cd ..    && $(MAKE) lib

lib:
	@echo ======================================== lib
	cd ../.. && $(MAKE) lib
	cd ..    && $(MAKE) lib


# depend on header, in case struct magma_tally3_opts changes
$(ALLOBJ): $(MAGMA_tally3_DIR)/testing/testings.h


$(EXE): $(LIBMAGMA_tally3) $(LIBMAGMA_tally3_SPARSE) $(LIBTEST)

clean: cleanobj cleanexe

cleanobj:
	rm -f *.$(o_ext) *~

cleanexe:
	rm -f $(EXE)

cleanall: clean cleanexe

testing_% : testing_%.$(o_ext)
	$(CXX) $(LDFLAGS) $(RPATH) $(NOMAIN) $< -o $@ \
	$(LIBTEST) -L$(MAGMA_tally3_DIR)/lib -lmagma_tally3_sparse -lmagma_tally3 \
	$(LIBDIR) \
	$(LIB)

run_% : run_%.$(o_ext)
	$(CXX) $(LDFLAGS) $(RPATH) $(NOMAIN) $< -o $@ \
	$(LIBTEST) -L$(MAGMA_tally3_DIR)/lib -lmagma_tally3_sparse -lmagma_tally3 \
	$(LIBDIR) \
	$(LIB)

# keep intermediate object files, to speed up re-compiling when library changes
.SECONDARY: $(ALLOBJ)
