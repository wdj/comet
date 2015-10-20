/*---------------------------------------------------------------------------*/
/*!
 * \file   compute_metrics_ccc.c
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Functions for computing CCC metrics.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include <stdlib.h>

#include "env.h"
#include "vectors.h"
#include "metrics.h"
#include "compute_metrics_ccc.h"

/*===========================================================================*/

void gm_compute_metrics_ccc_2way_cpu(GMMetrics* metrics,
                                     GMVectors* vectors,
                                     GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  if (env->all2all) {
  
    /*---Initializations---*/

    int data_type = gm_data_type_from_metric_type(env->metric_type, env);

    /*--Create size-2 circular buffer of vectors objects for send/recv---*/
    
    GMVectors vectors_0 = GMVectors_null();
    GMVectors vectors_1 = GMVectors_null();
    
    GMVectors_create(&vectors_0, data_type, vectors->num_field,
                     vectors->num_vector_local, env);
    GMVectors_create(&vectors_1, data_type, vectors->num_field,
                     vectors->num_vector_local, env);  

    /*---Result matrix is diagonal block and half the blocks to the right
         (including wraparound to left side of matrix when appropriate).
         For even number of procs, block rows of lower half of matrix
         have one less block to make correct count---*/

    const int num_step = 1 + (env->num_proc / 2);

    /*---Loop over steps of circular shift of vectors objects---*/

    int step_num = 0;
    for (step_num = 0; step_num < num_step; ++step_num) {
    
      const GMBool is_first_step = step_num == 0;
      const GMBool is_last_step = step_num != num_step - 1;

      /*---The proc that owns the "right-side" vectors for the minproduct---*/
    
      const int j_proc = (env->proc_num + step_num) % env->num_proc;

      /*---Left- and right-side vecs, also right-side vecs for next step---*/

      GMVectors* vectors_left = vectors;
      GMVectors* vectors_right = 
          is_first_step ? vectors : vectors_01[step_num % 2];
      GMVectors* vectors_right_next = vectors_01[(step_num + 1) % 2];

      /*---Initialize sends/recvs---*/
 
      MPI_Request mpi_requests[2];
 
      int mpi_code = 0;
      if (mpi_code) {
      }

      const int proc_up = (env->proc_num + 1) % env->num_proc;
      const int proc_dn = (env->proc_num -1 + env->num_proc) % env->num_proc;

      const int mpi_tag = 0;

      /*---Initiate sends/recvs for vecs needed on next step---*/

      if (is_last_step) {
        mpi_code = MPI_Isend(
            (void*)vectors_right->data, vectors_right->num_dataval_local,
            GM_MPI_FlOAT, proc_dn, mpi_tag, env->mpi_comm, &(mpi_requests[0]));
        GMAssert(mpi_code == MPI_SUCCESS);
      } /*---if step_num---*/

      /*---Compute frequency factors F_i = 1 - f_i/q---*/
      /*---We need to add a field in vectors that holds the frequencies and alleles---*/


  GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
}

/*===========================================================================*/

void gm_compute_metrics_ccc_2way_gpu(GMMetrics* metrics,
                                     GMVectors* vectors,
                                     GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  GMInsist(env, (!env->all2all) ? "Unimplemented." : 0);

  GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
}

/*===========================================================================*/

void gm_compute_metrics_ccc_3way_cpu(GMMetrics* metrics,
                                     GMVectors* vectors,
                                     GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  GMInsist(env, (!env->all2all) ? "Unimplemented." : 0);

  GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
}

/*===========================================================================*/

void gm_compute_metrics_ccc_3way_gpu(GMMetrics* metrics,
                                     GMVectors* vectors,
                                     GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  GMInsist(env, (!env->all2all) ? "Unimplemented." : 0);

  GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
}

/*===========================================================================*/

/*---------------------------------------------------------------------------*/
