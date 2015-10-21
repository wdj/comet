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
    
    GMFloat* freq_onproc = GMFloat_malloc(2 * metrics->num_vector_local);
    GMFloat* freq_offproc = GMFloat_malloc(2 * metrics->num_vector_local);
    GMFloat* allele_onproc = (char *) malloc(2 * metrics->num_vector_local * sizeof(char)); 
    GMFloat* allele_offproc = (char *) malloc(2 * metrics->num_vector_local * sizeof(char)); 

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
      
 
        mpi_code = MPI_Irecv((void*)vectors_right_next->data,
                             vectors_right_next->num_dataval_local, GM_MPI_FLOAT,
                             proc_up, mpi_tag, env->mpi_comm, &(mpi_requests[1]));
        GMAssert(mpi_code == MPI_SUCCESS);

      } /*---if step_num---*/

      /*---Compute metrics---*/
      
      const GMBool compute_triang_only = is_first_step;

      const GMBool skipped_last_block_lower_half = env->num_proc % 2 == 0 &&
                        ( 2 * env->proc_num >= env->num_proc ) && is_last_step;
  
      if ( ! skipped_last_block_lower_half ) {
        int j = 0;
        for (j = 0; j < metrics->num_vector_local; ++j) {
          const int i_max = compute_triang_only ? j : metrics->num_vector_local;
          int i = 0;
          for (i = 0; i < i_max; ++i) {
            /*---Add up number of individuals with each relationship---*/
            GMFloat* tally = GMFloat_malloc(16);
            int index1 = 0;
            int index2 = 0;
            for (k = 0; k < numInd; k++) {
              index1 = (int)GMVectors_float_get(vectors, field, i, env);
              index2 = (int)GMVectors_float_get(vectors, field, j, env);
              tally[(index1 + 4*index2)]++;
            }

            /*---Adjust proportionate contributions of each relationship---*/
            tally[5] /= 4.0; // both heterozygous
            tally[4] /= 2.0; // one heterozygous, the other homozygous
            tally[1] /= 2.0; // one heterozygous, the other homozygous
            tally[9] /= 2.0; // one heterozygous, the other homozygous
            tally[6] /= 2.0; // one heterozygous, the other homozygous

            /*---Count how many individuals have no missing data---*/
            int noMissing = 0;
            int col = 0;
            int row = 0;
            for (col = 0; col < 3; col++) {
              for (row = 0; row < 3; row++) {
                noMissing += (int)tally[row+4*col];
              }
            }
            
            /*---Initialize the four output metrics: ll, lh, hl, and hh--*/
            /*---TO DO: How do we incorporate these 4 output metrics into the metric struct?---*/
            int numSnps = metrics->num_vector_local;
            ll = (float *) malloc(numSnps*numSnps*sizeof(float));
  	    lh = (float *) malloc(numSnps*numSnps*sizeof(float));
            hl = (float *) malloc(numSnps*numSnps*sizeof(float));
            hh = (float *) malloc(numSnps*numSnps*sizeof(float));

            /*---Compute four relationship values---*/
            ll[j+numSnps*i] = tally[0] + tally[4] + tally[2] + tally[5];
            lh[j+numSnps*i] = tally[4] + tally[8] + tally[5] + tally[9];
            hl[j+numSnps*i] = tally[1] + tally[5] + tally[2] + tally[6];
            hh[j+numSnps*i] = tally[5] + tally[9] + tally[6] + tally[10];

            /*---Find average by dividing by number of individuals---*/
            ll[j+numSnps*i] /= (float)noMissing;
            lh[j+numSnps*i] /= (float)noMissing;
            hl[j+numSnps*i] /= (float)noMissing;
            hh[j+numSnps*i] /= (float)noMissing;
          }
        }
      }
    }/*---for step_num---*/ 

  } else /*---if (! env->all2all )---*/ {
    
    int i = 0;
    for (i = 0; i < metrics->num_vector_local; ++i) {
      int j = 0;
      for (j = i + 1; j < metrics->num_vector_local; ++j) {

        /*---Add up number of individuals with each relationship---*/
        GMFloat* tally = GMFloat_malloc(16);
        int index1 = 0;
        int index2 = 0;
        for (k = 0; k < numInd; k++) {
          index1 = (int)GMVectors_float_get(vectors, field, i, env);
          index2 = (int)GMVectors_float_get(vectors, field, j, env); 
          tally[(index1 + 4*index2)]++;
        }
         
        /*---Adjust proportionate contributions of each relationship---*/
        tally[5] /= 4.0; // both heterozygous
        tally[4] /= 2.0; // one heterozygous, the other homozygous
        tally[1] /= 2.0; // one heterozygous, the other homozygous
        tally[9] /= 2.0; // one heterozygous, the other homozygous
        tally[6] /= 2.0; // one heterozygous, the other homozygous

        /*---Count how many individuals have no missing data---*/
        int noMissing = 0;
        int col = 0;
        int row = 0;
        for (col = 0; col < 3; col++) {
          for (row = 0; row < 3; row++) {
            noMissing += (int)tally[row+4*col];
          }
        }

        /*---Initialize the four output metrics: ll, lh, hl, and hh--*
        /*---TO DO: How do we incorporate these 4 output metrics into the metric struct?---*/
        int numSnps = metrics->num_vector_local;
        ll = (float *) malloc(numSnps*numSnps*sizeof(float));
        lh = (float *) malloc(numSnps*numSnps*sizeof(float));
        hl = (float *) malloc(numSnps*numSnps*sizeof(float));
        hh = (float *) malloc(numSnps*numSnps*sizeof(float));

        /*---Compute four relationship values---*/
        ll[j+numSnps*i] = tally[0] + tally[4] + tally[2] + tally[5];
        lh[j+numSnps*i] = tally[4] + tally[8] + tally[5] + tally[9];
        hl[j+numSnps*i] = tally[1] + tally[5] + tally[2] + tally[6];
        hh[j+numSnps*i] = tally[5] + tally[9] + tally[6] + tally[10];

        /*---Find average by dividing by number of individuals---*/
        ll[j+numSnps*i] /= (float)noMissing;
        lh[j+numSnps*i] /= (float)noMissing;
        hl[j+numSnps*i] /= (float)noMissing;
        hh[j+numSnps*i] /= (float)noMissing;
 
        /*---Multiply by frequency factors---*/
        ll[j+numSnps*i] *= freq[i] * freq[j]; 
        lh[j+numSnps*i] *= freq[i] * freq[j+numSnps];
        hl[j+numSnps*i] *= freq[i+numSnps] * freq[j];
        hh[j+numSnps*i] *= freq[i+numSnps] * freq[j+numSnps];
      }
    }

  }/*---if (env->all2all )---*/

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
