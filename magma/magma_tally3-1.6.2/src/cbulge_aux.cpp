/*
 * Copyright (c) 2011      The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 *
 *     @author Azzam Haidar
 *     @author Stan Tomov
 *
 *     @generated from zbulge_aux.cpp normal z -> c, Fri Jan 30 19:00:18 2015
 *
 */

#include "common_magma_tally3.h"
#include "magma_tally3_timer.h"


//////////////////////////////////////////////////////////////
//          CSTEDC          Divide and Conquer for tridiag
//////////////////////////////////////////////////////////////
extern "C" void  magma_tally3_cstedc_withZ(magma_tally3_vec_t JOBZ, magma_tally3_int_t N, float *D, float * E, magma_tally3FloatComplex *Z, magma_tally3_int_t LDZ)
{
    magma_tally3FloatComplex *WORK;
    float *RWORK;
    magma_tally3_int_t *IWORK;
    magma_tally3_int_t LWORK, LIWORK, LRWORK;
    magma_tally3_int_t INFO;
    
    // use log() as log2() is not defined everywhere (e.g., Windows)
    const float log_2 = 0.6931471805599453;
    if (JOBZ == Magma_tally3Vec) {
        LWORK  = N*N;
        LRWORK = 1 + 3*N + 3*N*((magma_tally3_int_t)(log( (float)N )/log_2) + 1) + 4*N*N + 256*N;
        LIWORK = 6 + 6*N + 6*N*((magma_tally3_int_t)(log( (float)N )/log_2) + 1) + 256*N;
    } else if (JOBZ == Magma_tally3IVec) {
        LWORK  = N;
        LRWORK = 2*N*N + 4*N + 1 + 256*N;
        LIWORK = 256*N;
    } else if (JOBZ == Magma_tally3NoVec) {
        LWORK  = N;
        LRWORK = 256*N + 1;
        LIWORK = 256*N;
    } else {
        printf("ERROR JOBZ %c\n", JOBZ);
        exit(-1);
    }
    
    magma_tally3_smalloc_cpu( &RWORK, LRWORK );
    magma_tally3_cmalloc_cpu( &WORK,  LWORK  );
    magma_tally3_imalloc_cpu( &IWORK, LIWORK );
    
    lapackf77_cstedc( lapack_vec_const_tally3(JOBZ), &N, D, E, Z, &LDZ, WORK, &LWORK, RWORK, &LRWORK, IWORK, &LIWORK, &INFO);
    
    if (INFO != 0) {
        printf("=================================================\n");
        printf("CSTEDC ERROR OCCURED. HERE IS INFO %d \n ", (int) INFO);
        printf("=================================================\n");
        //assert(INFO == 0);
    }
    
    magma_tally3_free_cpu( IWORK );
    magma_tally3_free_cpu( WORK  );
    magma_tally3_free_cpu( RWORK );
}
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//          CSTEDC          Divide and Conquer for tridiag
//////////////////////////////////////////////////////////////
extern "C" void  magma_tally3_cstedx_withZ(magma_tally3_int_t N, magma_tally3_int_t NE, float *D, float * E, magma_tally3FloatComplex *Z, magma_tally3_int_t LDZ)
{
    float *RWORK;
    float *dwork;
    magma_tally3_int_t *IWORK;
    magma_tally3_int_t LIWORK, LRWORK;
    magma_tally3_int_t INFO;
    
    //LWORK  = N;
    LRWORK = 2*N*N + 4*N + 1 + 256*N;
    LIWORK = 256*N;
    
    magma_tally3_smalloc_cpu( &RWORK, LRWORK );
    magma_tally3_imalloc_cpu( &IWORK, LIWORK );
    
    if (MAGMA_tally3_SUCCESS != magma_tally3_smalloc( &dwork, 3*N*(N/2 + 1) )) {
        printf("=================================================\n");
        printf("CSTEDC ERROR OCCURED IN CUDAMALLOC\n");
        printf("=================================================\n");
        return;
    }
    printf("using magma_tally3_cstedx\n");

    magma_tally3_timer_t time=0;
    timer_start( time );

    magma_tally3_range_t job = Magma_tally3RangeI;
    if (NE == N)
        job = Magma_tally3RangeAll;

    magma_tally3_cstedx(job, N, 0., 0., 1, NE, D, E, Z, LDZ, RWORK, LRWORK, IWORK, LIWORK, dwork, &INFO);

    if (INFO != 0) {
        printf("=================================================\n");
        printf("CSTEDC ERROR OCCURED. HERE IS INFO %d \n ", (int) INFO);
        printf("=================================================\n");
        //assert(INFO == 0);
    }

    timer_stop( time );
    timer_printf( "time zstevx = %6.2f\n", time );

    magma_tally3_free( dwork );
    magma_tally3_free_cpu( IWORK );
    magma_tally3_free_cpu( RWORK );
}
