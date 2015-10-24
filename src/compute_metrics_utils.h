/*---------------------------------------------------------------------------*/
/*!
 * \file   compute_metrics_utils.h
 * \author Wayne Joubert
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Functions for computing metrics, utilities, header.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _compute_metrics_utils_h_
#define _compute_metrics_utils_h_

#include "env.h"
#include "vector_sums.h"
#include "vectors.h"
#include "metrics.h"

/*===========================================================================*/

void gm_compute_float_vector_sums(GMVectors* vectors,
                                  GMFloat* __restrict__ vector_sums,
                                  GMEnv* env);

void gm_compute_vector_sums(GMVectors* vectors,
                            GMVectorSums* vector_sums,
                            GMEnv* env);

/*----------*/

void gm_magma_initialize(GMEnv* env);

void gm_magma_finalize(GMEnv* env);

/*----------*/

GMMirroredPointer gm_malloc_magma(size_t n, GMEnv* env);

void gm_free_magma(GMMirroredPointer* p, GMEnv* env);

/*----------*/

MPI_Request gm_send_vectors_start(GMVectors* vectors, 
                                  int proc_num,
                                  GMEnv* env);

MPI_Request gm_recv_vectors_start(GMVectors* vectors, 
                                  int proc_num,
                                  GMEnv* env);

void gm_send_vectors_wait(MPI_Request* mpi_request, GMEnv* env);

void gm_recv_vectors_wait(MPI_Request* mpi_request, GMEnv* env);

/*----------*/

void gm_set_vectors_start(GMVectors* vectors,
                          GMMirroredPointer* vectors_buf,
                          GMEnv* env);

void gm_set_vectors_wait(GMEnv* env);

void gm_get_metrics_start(GMMetrics* metrics,
                          GMMirroredPointer* metrics_buf,
                          GMEnv* env);

void gm_get_metrics_wait(GMEnv* env);

/*----------*/

void gm_vectors_to_buf(GMVectors* vectors,
                       GMMirroredPointer* vectors_buf,
                       GMEnv* env);


/*----------*/

void gm_compute_numerators_start( GMVectors* vectors_left,
                                  GMVectors* vectors_right,
                                  GMMetrics* numerators,
                                  GMMirroredPointer* vectors_left_buf,
                                  GMMirroredPointer* vectors_right_buf,
                                  GMMirroredPointer* numerators_buf,
                                  int j_proc,
                                  GMBool compute_triang_only,
                                  GMEnv* env);

void gm_compute_numerators_wait(GMEnv* env);

void gm_compute_czekanowski_combine(GMMetrics* metrics,
                                    GMMirroredPointer* metrics_buf,
                                    GMFloat* __restrict__ vector_sums_left,
                                    GMFloat* __restrict__ vector_sums_right,
                                    int j_proc,
                                    GMBool compute_triang_only,
                                    GMEnv* env);

void gm_compute_combine(GMMetrics* metrics,
                        GMMirroredPointer* metrics_buf,
                        GMVectorSums* vector_sums_left,
                        GMVectorSums* vector_sums_right,
                        int j_proc, 
                        GMBool do_compute_triang_only,
                        GMEnv* env);

/*===========================================================================*/

#endif /*---_compute_metrics_utils_h---*/

/*---------------------------------------------------------------------------*/
