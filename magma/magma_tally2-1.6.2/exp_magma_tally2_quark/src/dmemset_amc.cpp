/* 
    -- MAGMA_tally2 (version 1.6.1) -- 
       Univ. of Tennessee, Knoxville 
       Univ. of California, Berkeley 
       Univ. of Colorado, Denver 
       Sept 2013 
 
       @author: Simplice Donfack 
*/

#include "common_magma_tally2.h"

#include "magma_tally2_amc.h"

#include "schedule.h"

/*Initialize the matrix in parallel. This may prevent to use numactl --interleave*/
void magma_tally2_dmemset_amc(double *ptr, double value, int n,  magma_tally2_int_t chunck, int P){

    int i, nsize, info;

    /* Check arguments */ 
    info = 0; 
    if (n < 0) 
        info = -1; 
    else if (chunck < 0) 
        info = -2; 
 
    if (info != 0) { 
        magma_tally2_xerbla( __func__, -info ); 
        return ; 
    } 
 
    /* Quick return if possible */ 
    if (n == 0 || chunck == 0) 
        return;
    
    /*Initialize the scheduler*/ 
#ifdef NGPU
    magma_tally2_schedule_init(P, NGPU);
#else
    magma_tally2_schedule_init(P, 1);
#endif
     

    for(i=0;i<n;i+=chunck){
        nsize = min(n-i, chunck);
        magma_tally2_insert_core_dmemset(&ptr[i], value, nsize);
    }

    /*Wait for all thread termination*/
    magma_tally2_schedule_barrier();

    /*Shutdown the scheduler*/
     magma_tally2_schedule_delete();
}

