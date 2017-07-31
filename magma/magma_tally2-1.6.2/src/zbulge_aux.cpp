/*
 * Copyright (c) 2011      The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 *
 *     @author Azzam Haidar
 *     @author Stan Tomov
 *
 *     @precisions normal z -> c
 *
 */

#include "common_magma_tally2.h"
#include "magma_tally2_timer.h"


//////////////////////////////////////////////////////////////
//          ZSTEDC          Divide and Conquer for tridiag
//////////////////////////////////////////////////////////////
extern "C" void  magma_tally2_zstedc_withZ(magma_tally2_vec_t JOBZ, magma_tally2_int_t N, double *D, double * E, magma_tally2DoubleComplex *Z, magma_tally2_int_t LDZ)
{
    magma_tally2DoubleComplex *WORK;
    double *RWORK;
    magma_tally2_int_t *IWORK;
    magma_tally2_int_t LWORK, LIWORK, LRWORK;
    magma_tally2_int_t INFO;
    
    // use log() as log2() is not defined everywhere (e.g., Windows)
    const double log_2 = 0.6931471805599453;
    if (JOBZ == Magma_tally2Vec) {
        LWORK  = N*N;
        LRWORK = 1 + 3*N + 3*N*((magma_tally2_int_t)(log( (double)N )/log_2) + 1) + 4*N*N + 256*N;
        LIWORK = 6 + 6*N + 6*N*((magma_tally2_int_t)(log( (double)N )/log_2) + 1) + 256*N;
    } else if (JOBZ == Magma_tally2IVec) {
        LWORK  = N;
        LRWORK = 2*N*N + 4*N + 1 + 256*N;
        LIWORK = 256*N;
    } else if (JOBZ == Magma_tally2NoVec) {
        LWORK  = N;
        LRWORK = 256*N + 1;
        LIWORK = 256*N;
    } else {
        printf("ERROR JOBZ %c\n", JOBZ);
        exit(-1);
    }
    
    magma_tally2_dmalloc_cpu( &RWORK, LRWORK );
    magma_tally2_zmalloc_cpu( &WORK,  LWORK  );
    magma_tally2_imalloc_cpu( &IWORK, LIWORK );
    
    lapackf77_zstedc( lapack_vec_const_tally2(JOBZ), &N, D, E, Z, &LDZ, WORK, &LWORK, RWORK, &LRWORK, IWORK, &LIWORK, &INFO);
    
    if (INFO != 0) {
        printf("=================================================\n");
        printf("ZSTEDC ERROR OCCURED. HERE IS INFO %d \n ", (int) INFO);
        printf("=================================================\n");
        //assert(INFO == 0);
    }
    
    magma_tally2_free_cpu( IWORK );
    magma_tally2_free_cpu( WORK  );
    magma_tally2_free_cpu( RWORK );
}
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//          ZSTEDC          Divide and Conquer for tridiag
//////////////////////////////////////////////////////////////
extern "C" void  magma_tally2_zstedx_withZ(magma_tally2_int_t N, magma_tally2_int_t NE, double *D, double * E, magma_tally2DoubleComplex *Z, magma_tally2_int_t LDZ)
{
    double *RWORK;
    double *dwork;
    magma_tally2_int_t *IWORK;
    magma_tally2_int_t LIWORK, LRWORK;
    magma_tally2_int_t INFO;
    
    //LWORK  = N;
    LRWORK = 2*N*N + 4*N + 1 + 256*N;
    LIWORK = 256*N;
    
    magma_tally2_dmalloc_cpu( &RWORK, LRWORK );
    magma_tally2_imalloc_cpu( &IWORK, LIWORK );
    
    if (MAGMA_tally2_SUCCESS != magma_tally2_dmalloc( &dwork, 3*N*(N/2 + 1) )) {
        printf("=================================================\n");
        printf("ZSTEDC ERROR OCCURED IN CUDAMALLOC\n");
        printf("=================================================\n");
        return;
    }
    printf("using magma_tally2_zstedx\n");

    magma_tally2_timer_t time=0;
    timer_start( time );

    magma_tally2_range_t job = Magma_tally2RangeI;
    if (NE == N)
        job = Magma_tally2RangeAll;

    magma_tally2_zstedx(job, N, 0., 0., 1, NE, D, E, Z, LDZ, RWORK, LRWORK, IWORK, LIWORK, dwork, &INFO);

    if (INFO != 0) {
        printf("=================================================\n");
        printf("ZSTEDC ERROR OCCURED. HERE IS INFO %d \n ", (int) INFO);
        printf("=================================================\n");
        //assert(INFO == 0);
    }

    timer_stop( time );
    timer_printf( "time zstevx = %6.2f\n", time );

    magma_tally2_free( dwork );
    magma_tally2_free_cpu( IWORK );
    magma_tally2_free_cpu( RWORK );
}
