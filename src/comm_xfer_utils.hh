//-----------------------------------------------------------------------------
/*!
 * \file   comm_xfer_utils.hh
 * \author Wayne Joubert
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Utilities for communication and CPU/GPU transfer of vectors, metrics.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _gm_comm_xfer_utils_hh_
#define _gm_comm_xfer_utils_hh_

#include "env.hh"
#include "mirrored_buf.hh"
#include "vectors.hh"
#include "metrics.hh"

//=============================================================================

MPI_Request gm_send_vectors_start(GMVectors* vectors, int proc_num, int mpi_tag,
                                  GMEnv* env);

MPI_Request gm_recv_vectors_start(GMVectors* vectors, int proc_num, int mpi_tag,
                                  GMEnv* env);

void gm_send_vectors_wait(MPI_Request* mpi_request, GMEnv* env);

void gm_recv_vectors_wait(MPI_Request* mpi_request, GMEnv* env);

//--------------------

void gm_reduce_metrics(GMMetrics* metrics,
                       GMMirroredBuf* metrics_buf_target,
                       GMMirroredBuf* metrics_buf_source,
                       GMEnv* env);

MPI_Request gm_reduce_metrics_start(GMMetrics* metrics,
                                    GMMirroredBuf* metrics_buf_target,
                                    GMMirroredBuf* metrics_buf_source,
                                    GMEnv* env);

void gm_reduce_metrics_wait(MPI_Request* mpi_request, GMEnv* env);

//--------------------

void gm_set_vectors_start(GMVectors* vectors,
                          GMMirroredBuf* vectors_buf,
                          GMEnv* env);

void gm_set_vectors_wait(GMEnv* env);

void gm_get_metrics_start(GMMetrics* metrics,
                          GMMirroredBuf* metrics_buf,
                          GMEnv* env);

void gm_get_metrics_wait(GMMetrics* metrics,
                         GMMirroredBuf* metrics_buf,
                         GMEnv* env);

//=============================================================================

#endif // _gm_comm_xfer_utils_hh_

//-----------------------------------------------------------------------------
