#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA_minproduct (version 1.6.1) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      @date January 2015
#//////////////////////////////////////////////////////////////////////////////

MAGMA_minproduct_DIR = .
include ./Makefile.internal

all: lib test

lib: libmagma_minproduct libquark

libmagma_minproduct:
	( cd control     && $(MAKE) )
	( cd src         && $(MAKE) )

libquark:
	( cd quark       && $(MAKE) )

test:
	cp ../lib/libmagma_minproductblas.a lib/libmagma_minproductblas.a
	( cd testing     && $(MAKE) )

clean:
	( cd include     && $(MAKE) clean )
	( cd control     && $(MAKE) clean )
	( cd src         && $(MAKE) clean )
	( cd testing     && $(MAKE) clean )
	( cd quark       && $(MAKE) clean )

cleanall:
	( cd include     && $(MAKE) cleanall )
	( cd control     && $(MAKE) cleanall )
	( cd src         && $(MAKE) cleanall )
	( cd testing     && $(MAKE) cleanall )
	( cd quark       && $(MAKE) cleanall )
	( cd lib         && rm -f *.a )
	$(MAKE) cleanall2

# cleanall2 is a dummy rule to run cleangen at the *end* of make cleanall, so
# .Makefile.gen files aren't deleted and immediately re-created. see Makefile.gen
cleanall2:
	@echo

dir:
	mkdir -p $(prefix)
	mkdir -p $(prefix)/include
	mkdir -p $(prefix)/lib
	mkdir -p $(prefix)/lib/pkgconfig

install: lib dir
#       MAGMA_minproduct
	cp $(MAGMA_minproduct_DIR)/include/*.h  $(prefix)/include
	cp $(LIBMAGMA_minproduct)               $(prefix)/lib
	cp $(LIBMAGMA_minproductBLAS)           $(prefix)/lib
#       QUARK
	cp $(QUARKDIR)/include/quark.h             $(prefix)/include
	cp $(QUARKDIR)/include/quark_unpack_args.h $(prefix)/include
	cp $(QUARKDIR)/include/icl_hash.h          $(prefix)/include
	cp $(QUARKDIR)/include/icl_list.h          $(prefix)/include
	cp $(QUARKDIR)/lib/libquark.a              $(prefix)/lib
#       pkgconfig
	cat $(MAGMA_minproduct_DIR)/lib/pkgconfig/magma_minproduct.pc | \
	    sed -e s:\__PREFIX:"$(prefix)":     | \
	    sed -e s:\__LIBEXT:"$(LIBEXT)":       \
	    > $(prefix)/lib/pkgconfig/magma_minproduct.pc
