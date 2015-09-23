/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Raffaele Solca

*/
#ifndef MAGMA_minproduct_AFFINITY_H
#define MAGMA_minproduct_AFFINITY_H

#ifdef MAGMA_minproduct_SETAFFINITY

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <sched.h>

#if __GLIBC_PREREQ(2,3)

class affinity_set
{
public:

    affinity_set();

    affinity_set(int cpu_nr);

    void add(int cpu_nr);

    int get_affinity();

    int set_affinity();

    void print_affinity(int id, const char* s);

    void print_set(int id, const char* s);

private:

    cpu_set_t set;
};

#else
#error "Affinity requires Linux glibc version >= 2.3.3, which isn't available. Remove -DMAGMA_minproduct_SETAFFINITY from CFLAGS in make.inc."
#endif

#endif  // MAGMA_minproduct_SETAFFINITY

#endif  // MAGMA_minproduct_AFFINITY_H
