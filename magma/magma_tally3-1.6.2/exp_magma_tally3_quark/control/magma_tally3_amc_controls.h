/* 
    -- MAGMA_tally3 (version 1.6.1) -- 
       Univ. of Tennessee, Knoxville 
       Univ. of California, Berkeley 
       Univ. of Colorado, Denver 
       May 2013 
 
       @author: Simplice Donfack 
 
*/

/*The routines control the initialisation of the routines called in magma_tally3_amc*/
#ifndef MAGMA_tally3_AMC_CONTROLS_H
#define MAGMA_tally3_AMC_CONTROLS_H

int magma_tally3_amc_init(int P, double dcpu, int Pr, int nb);
int magma_tally3_amc_finalize();
#endif

