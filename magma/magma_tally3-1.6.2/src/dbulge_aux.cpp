/*
 * Copyright (c) 2011      The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 *
 *     @author Azzam Haidar
 *     @author Stan Tomov
 *
 *     @precisions normal d -> s
 *
 */

#include "common_magma_tally3.h"
#include "magma_tally3_timer.h"
#include "magma_tally3_dbulgeinc.h"

//////////////////////////////////////////////////////////////
//          DSTEDC          Divide and Conquer for tridiag
//////////////////////////////////////////////////////////////
extern "C" void  magma_tally3_dstedc_withZ(magma_tally3_vec_t JOBZ, magma_tally3_int_t N, double *D, double * E, double *Z, magma_tally3_int_t LDZ)
{
    double *WORK;
    magma_tally3_int_t *IWORK;
    magma_tally3_int_t LWORK, LIWORK;
    magma_tally3_int_t INFO;
    
    // use log() as log2() is not defined everywhere (e.g., Windows)
    const double log_2 = 0.6931471805599453;
    if (JOBZ == Magma_tally3Vec) {
        LWORK  = 1 + 3*N + 3*N*((magma_tally3_int_t)(log( (double)N )/log_2) + 1) + 4*N*N + 256*N;
        LIWORK = 6 + 6*N + 6*N*((magma_tally3_int_t)(log( (double)N )/log_2) + 1) + 256*N;
    } else if (JOBZ == Magma_tally3IVec) {
        LWORK  = 2*N*N + 256*N + 1;
        LIWORK = 256*N;
    } else if (JOBZ == Magma_tally3NoVec) {
        LWORK  = 256*N + 1;
        LIWORK = 256*N;
    } else {
        printf("ERROR JOBZ %c\n", JOBZ);
        exit(-1);
    }
    
    magma_tally3_dmalloc_cpu( &WORK,  LWORK  );
    magma_tally3_imalloc_cpu( &IWORK, LIWORK );
    
    lapackf77_dstedc( lapack_vec_const_tally3(JOBZ), &N, D, E, Z, &LDZ, WORK, &LWORK, IWORK, &LIWORK, &INFO);
    
    if (INFO != 0) {
        printf("=================================================\n");
        printf("DSTEDC ERROR OCCURED. HERE IS INFO %d \n ", (int) INFO);
        printf("=================================================\n");
        //assert(INFO == 0);
    }
    
    magma_tally3_free_cpu( IWORK );
    magma_tally3_free_cpu( WORK  );
}
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//          DSTEDX          Divide and Conquer for tridiag
//////////////////////////////////////////////////////////////
extern "C" void  magma_tally3_dstedx_withZ(magma_tally3_int_t N, magma_tally3_int_t NE, double *D, double * E, double *Z, magma_tally3_int_t LDZ)
{
    double *WORK;
    double *dwork;
    magma_tally3_int_t *IWORK;
    magma_tally3_int_t LWORK, LIWORK;
    magma_tally3_int_t INFO;
    
    LWORK  = N*N+4*N+1;
    LIWORK = 3 + 5*N;
    
    magma_tally3_dmalloc_cpu( &WORK,  LWORK  );
    magma_tally3_imalloc_cpu( &IWORK, LIWORK );
    
    if (MAGMA_tally3_SUCCESS != magma_tally3_dmalloc( &dwork, 3*N*(N/2 + 1) )) {
        printf("=================================================\n");
        printf("DSTEDC ERROR OCCURED IN CUDAMALLOC\n");
        printf("=================================================\n");
        return;
    }
    printf("using magma_tally3_dstedx\n");

    magma_tally3_timer_t time=0;
    timer_start( time );

    //magma_tally3_range_t job = Magma_tally3RangeI;
    //if (NE == N)
    //    job = Magma_tally3RangeAll;
    
    magma_tally3_dstedx(Magma_tally3RangeI, N, 0., 0., 1, NE, D, E, Z, LDZ, WORK, LWORK, IWORK, LIWORK, dwork, &INFO);
    
    if (INFO != 0) {
        printf("=================================================\n");
        printf("DSTEDC ERROR OCCURED. HERE IS INFO %d \n ", (int) INFO);
        printf("=================================================\n");
        //assert(INFO == 0);
    }

    timer_stop( time );
    timer_printf( "time dstedx = %6.2f\n", time );

    magma_tally3_free( dwork );
    magma_tally3_free_cpu( IWORK );
    magma_tally3_free_cpu( WORK  );
}
//////////////////////////////////////////////////////////////
