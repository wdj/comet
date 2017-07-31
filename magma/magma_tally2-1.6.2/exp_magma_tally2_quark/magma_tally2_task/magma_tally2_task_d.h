/* 
    -- MAGMA_tally2 (version 1.6.1) -- 
       Univ. of Tennessee, Knoxville 
       Univ. of California, Berkeley 
       Univ. of Colorado, Denver 
       May 2013 
 
       @author: Simplice Donfack 
*/
#ifndef MAGMA_tally2_TASK_D_H
#define MAGMA_tally2_TASK_D_H

#include "schedule.h"

void magma_tally2_task_dmalloc_pinned(Schedule* sched_obj);
void magma_tally2_task_dfree_pinned(Schedule* sched_obj);
void magma_tally2_task_dfree_pinned_index(Schedule* sched_obj);
void magma_tally2_task_dsetmatrix(Schedule* sched_obj);
void magma_tally2_task_dgetmatrix(Schedule* sched_obj);
void magma_tally2_task_dsetmatrix_transpose(Schedule* sched_obj);
void magma_tally2_task_dgetmatrix_transpose(Schedule* sched_obj);               
void magma_tally2_task_dlaswp(Schedule* sched_obj);
void magma_tally2_task_dtrsm(Schedule* sched_obj);
void magma_tally2_task_dgemm(Schedule* sched_obj);
void magma_tally2_task_update(Schedule* sched_obj);
#endif

