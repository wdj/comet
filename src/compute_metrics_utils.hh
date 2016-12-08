/*---------------------------------------------------------------------------*/
/*!
 * \file   compute_metrics_utils.hh
 * \author Wayne Joubert, James Nance
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Functions for computing metrics, utilities, header.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _gm_compute_metrics_utils_hh_
#define _gm_compute_metrics_utils_hh_

#include "env.hh"
#include "vector_sums.hh"
#include "vectors.hh"
#include "metrics.hh"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/

MPI_Request gm_send_vectors_start(GMVectors* vectors, int proc_num, int mpi_tag,
                                  GMEnv* env);

MPI_Request gm_recv_vectors_start(GMVectors* vectors, int proc_num, int mpi_tag,
                                  GMEnv* env);

void gm_send_vectors_wait(MPI_Request* mpi_request, GMEnv* env);

void gm_recv_vectors_wait(MPI_Request* mpi_request, GMEnv* env);

/*----------*/

void gm_allreduce_metrics(GMMetrics* metrics,
                          GMMirroredPointer* metrics_buf_target,
                          GMMirroredPointer* metrics_buf_source,
                          GMEnv* env);

/*----------*/

void gm_set_vectors_start(GMVectors* vectors,
                          GMMirroredPointer* vectors_buf,
                          GMEnv* env);

void gm_set_vectors_wait(GMEnv* env);

void gm_get_metrics_start(GMMetrics* metrics,
                          GMMirroredPointer* metrics_buf,
                          GMEnv* env);

void gm_get_metrics_wait(GMMetrics* metrics,
                         GMMirroredPointer* metrics_buf,
                         GMEnv* env);

/*----------*/

void gm_vectors_to_buf(GMMirroredPointer* vectors_buf,
                       GMVectors* vectors,
                       GMEnv* env);

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*===========================================================================*/

#endif /*---_gm_compute_metrics_utils_hh_---*/

/*---------------------------------------------------------------------------*/
