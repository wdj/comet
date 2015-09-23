/* 
    -- MAGMA_minproduct (version 1.6.1) -- 
       Univ. of Tennessee, Knoxville 
       Univ. of California, Berkeley 
       Univ. of Colorado, Denver 
       May 2013 
 
       @author: Simplice Donfack 
 
*/

/*The routines control the initialisation of the routines called in magma_minproduct_amc*/
#ifndef MAGMA_minproduct_AMC_CONTROLS_H
#define MAGMA_minproduct_AMC_CONTROLS_H

int magma_minproduct_amc_init(int P, double dcpu, int Pr, int nb);
int magma_minproduct_amc_finalize();
#endif

