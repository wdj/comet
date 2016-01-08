/* 
    -- MAGMA_tally4 (version 1.6.1) -- 
       Univ. of Tennessee, Knoxville 
       Univ. of California, Berkeley 
       Univ. of Colorado, Denver 
       May 2013 
 
       @author: Simplice Donfack 
*/
/*
       Asynchronous and optimized multicore capabilities in MAGMA_tally4 (AMC).
*/
#ifndef MAGMA_tally4_AMC_H
#define MAGMA_tally4_AMC_H



#if (dbglevel >=1)
#include "ca_dbg_tools.h" /*Enable tracing tools and debugging tools*/
#endif

/* Some math and utilities functions */
#include "magma_tally4_amc_utils.h"

/* Context and arguments */
#include "magma_tally4_amc_args.h"

/* Initialisation and controls */
#include "magma_tally4_amc_controls.h"

/* Scheduling takes place here */
#include "schedule.h"

/* CPU tasks insertion */
#include "magma_tally4_insert_core_d.h"
/* 1 GPU tasks insertion */
#include "magma_tally4_insert_d.h"
/* Multi tasks insertion */
#include "magma_tally4_insert_dev_d.h"

/* Fill block of memory in async_dmemset.cpp*/
void magma_tally4_dmemset_amc(double *ptr, double value, int n, magma_tally4_int_t chunck, int P);

/*LU factorization in async_dgetrf_rec_gpu.cpp*/

extern "C" magma_tally4_int_t magma_tally4_dgetrf_gpu_amc(
magma_tally4_int_t m, magma_tally4_int_t n, 
double *dA, magma_tally4_int_t dA_LD,
magma_tally4_int_t *ipiv, magma_tally4_int_t *info
);

extern "C" magma_tally4_int_t magma_tally4_dgetrf_gpu_work_amc(
magma_tally4_int_t m, magma_tally4_int_t n,  
double *dA, magma_tally4_int_t dA_LD, 
magma_tally4_int_t *ipiv, magma_tally4_int_t *info,
/*Workspace on the cpu side*/
double *AWORK, magma_tally4_int_t AWORK_LD, magma_tally4_int_t AWORK_n
);

extern "C" magma_tally4_int_t magma_tally4_dgetrf_gpu_amc_v2(
magma_tally4_int_t m, magma_tally4_int_t n, 
double *dA, magma_tally4_int_t dA_LD,
magma_tally4_int_t *ipiv, magma_tally4_int_t *info
);

extern "C" magma_tally4_int_t magma_tally4_dgetrf_gpu_work_amc_v2(
magma_tally4_int_t m, magma_tally4_int_t n,  
double *dA, magma_tally4_int_t dA_LD, 
magma_tally4_int_t *ipiv, magma_tally4_int_t *info,
/*Workspace on the cpu side*/
double *AWORK, magma_tally4_int_t AWORK_LD, magma_tally4_int_t AWORK_n
);


extern "C" magma_tally4_int_t magma_tally4_dgetrf_mgpu_amc(magma_tally4_int_t num_gpus, magma_tally4_int_t m, magma_tally4_int_t n,  
                 double **dlA, magma_tally4_int_t dlA_LD, 
                 magma_tally4_int_t *ipiv, magma_tally4_int_t *info);

extern "C" magma_tally4_int_t 
magma_tally4_dgetrf_mgpu_work_amc(magma_tally4_int_t num_gpus,
magma_tally4_int_t m, magma_tally4_int_t n,  
double **dlA, magma_tally4_int_t dlA_LD, 
magma_tally4_int_t *ipiv, magma_tally4_int_t *info,
/*workspace on the cpu side*/
double *AWORK, magma_tally4_int_t AWORK_LD, magma_tally4_int_t AWORK_n,
/*workspace on the gpu side*/
double **dlpanelT, magma_tally4_int_t dlpanelT_m, magma_tally4_int_t dlpanelT_n
); 

extern "C" magma_tally4_int_t magma_tally4_dgetrf_mgpu_amc_v2(magma_tally4_int_t num_gpus, magma_tally4_int_t m, magma_tally4_int_t n,  
                 double **dlA, magma_tally4_int_t dlA_LD, 
                 magma_tally4_int_t *ipiv, magma_tally4_int_t *info);

extern "C" magma_tally4_int_t 
magma_tally4_dgetrf_mgpu_work_amc_v2(magma_tally4_int_t num_gpus,
magma_tally4_int_t m, magma_tally4_int_t n,  
double **dlA, magma_tally4_int_t dlA_LD, 
magma_tally4_int_t *ipiv, magma_tally4_int_t *info,
/*workspace on the cpu side*/
double *AWORK, magma_tally4_int_t AWORK_LD, magma_tally4_int_t AWORK_n,
/*workspace on the gpu side*/
double **dlpanelT, magma_tally4_int_t dlpanelT_m, magma_tally4_int_t dlpanelT_n
); 



extern "C" magma_tally4_int_t magma_tally4_dgetrf_mgpu_amc_v3(magma_tally4_int_t num_gpus, magma_tally4_int_t m, magma_tally4_int_t n,  
                 double **dlA, magma_tally4_int_t dlA_LD, 
                 magma_tally4_int_t *ipiv, magma_tally4_int_t *info);

extern "C" magma_tally4_int_t magma_tally4_dgetrf_mgpu_work_amc_v3(magma_tally4_int_t num_gpus,
magma_tally4_int_t m, magma_tally4_int_t n,  
double **dlA, magma_tally4_int_t dlA_LD, 
magma_tally4_int_t *ipiv, magma_tally4_int_t *info,
/*workspace on the cpu side*/
double *AWORK, magma_tally4_int_t AWORK_LD, magma_tally4_int_t AWORK_n
);

#endif

