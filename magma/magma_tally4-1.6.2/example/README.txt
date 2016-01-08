Getting started with MAGMA_tally4.

This is a simple, standalone example to show how to use MAGMA_tally4, once it is
compiled. More involved examples for individual routines are in the testing
directory. The testing code include some extra utilities that we use for
testing, such as testings.h and libtest.a, which are not required to use MAGMA_tally4,
though you may use them if desired.

If you want CUBLAS functions, include that header first, before magma_tally4.h:

    #include "cublas_v2.h"

Include the MAGMA_tally4 header:

    #include "magma_tally4.h"

You may also need BLAS and LAPACK functions, which you can get with:

    #include "magma_tally4_lapack.h"

You can also use headers that came with your BLAS and LAPACK library, such as
Intel MKL. However, their definitions, while compatible, may not exactly match
ours. Especially their definition of the COMPLEX type will be different. We use
magma_tally4DoubleComplex, which is a typedef of cuDoubleComplex. You may need to cast
back-and-forth between definitions.

When MAGMA_tally4 was compiled, one of ADD_, NOCHANGE, or UPCASE was defined in
make.inc for how Fortran functions are name-mangled on your system. The most
common is ADD_, where Fortran adds an underscore after the function. Usually,
you can tell the convention by using nm to examine your BLAS library:

    bunsen ~> nm /mnt/scratch/openblas/lib/libopenblas.a | grep -i dnrm2
    dnrm2.o:
    0000000000000000 T dnrm2_
                     U dnrm2_k

Since dnrm2_ has an underscore, use ADD_. Then add the include paths. For
example, to compile your .c file:

    gcc -DADD_ -DHAVE_CUBLAS \
        -I$(MAGMA_tally4DIR)/include -I$(CUDADIR)/include \
        -c -o example.c

where $(MAGMA_tally4DIR) and $(CUDADIR) are set to where MAGMA_tally4 and CUDA are installed,
respectively.

To link, add the library paths and necessary libraries. For example with OpenBLAS:

    gcc -L$(MAGMA_tally4DIR)/lib -L$(CUDADIR)/lib -L$(OPENBLASDIR) \
        -lmagma_tally4 -lcublas -lcudart -lopenblas \
        -o example example.o

The Makefile provided in this directory is a starting point for compiling. You
will need to adjust the MAGMA_tally4_CFLAGS and MAGMA_tally4_LIBS to reflect your system. See
the MAGMA_tally4 make.inc file for which libraries are required for BLAS & LAPACK.

Alternatively, you can use pkg-config to get MAGMA_tally4's CFLAGS and LIBS
automatically. To use pkg-config, install MAGMA_tally4 with 'make install', add MAGMA_tally4
to $PKG_CONFIG_PATH, e.g. with csh,

    setenv PKG_CONFIG_PATH ${PKG_CONFIG_PATH}:/usr/local/magma_tally4/lib/pkgconfig

then use 'pkg-config --cflags magma_tally4' and 'pkg-config --libs magma_tally4' as shown in
the Makefile.
