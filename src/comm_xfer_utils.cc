//-----------------------------------------------------------------------------
/*!
 * \file   comm_xfer_utils.cc
 * \author Wayne Joubert
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Communication, host/device transfer utilities.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#include "env.hh"
#include "mirrored_buf.hh"
#include "linalg.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "comm_xfer_utils.hh"

//=============================================================================
// Start/end MPI send/receive of vectors data

MPI_Request gm_send_vectors_start(GMVectors* vectors,
                                  int proc_num,
                                  int mpi_tag,
                                  GMEnv* env) {
  GMInsist(vectors && env);
  GMInsist(proc_num >= 0 && proc_num < GMEnv_num_proc_repl_vector(env));

  MPI_Request mpi_request;

  const MPI_Datatype mpi_type = gm_mpi_type(env);

  const int mpi_code =
      MPI_Isend((void*)vectors->data, vectors->num_packedval_local, mpi_type,
                proc_num, mpi_tag, GMEnv_mpi_comm_repl_vector(env),
                &mpi_request);
  GMInsist(mpi_code == MPI_SUCCESS);

  return mpi_request;
}

//-----------------------------------------------------------------------------

MPI_Request gm_recv_vectors_start(GMVectors* vectors,
                                  int proc_num,
                                  int mpi_tag,
                                  GMEnv* env) {
  GMInsist(vectors && env);
  GMInsist(proc_num >= 0 && proc_num < GMEnv_num_proc_repl_vector(env));

  MPI_Request mpi_request;

  const MPI_Datatype mpi_type = gm_mpi_type(env);

  const int mpi_code =
      MPI_Irecv((void*)vectors->data, vectors->num_packedval_local, mpi_type,
                proc_num, mpi_tag, GMEnv_mpi_comm_repl_vector(env),
                &mpi_request);
  GMInsist(mpi_code == MPI_SUCCESS);

  return mpi_request;
}

//-----------------------------------------------------------------------------

void gm_send_vectors_wait(MPI_Request* mpi_request, GMEnv* env) {
  GMInsist(mpi_request && env);

  MPI_Status mpi_status;

  const int mpi_code = MPI_Wait(mpi_request, &mpi_status);
  GMInsist(mpi_code == MPI_SUCCESS);
}

//-----------------------------------------------------------------------------

void gm_recv_vectors_wait(MPI_Request* mpi_request, GMEnv* env) {
  GMInsist(mpi_request && env);

  MPI_Status mpi_status;

  const int mpi_code = MPI_Wait(mpi_request, &mpi_status);
  GMInsist(mpi_code == MPI_SUCCESS);
}

//=============================================================================
// MPI reduce operations

void gm_reduce_metrics(GMMetrics* metrics,
                       GMMirroredBuf* metrics_buf_target,
                       GMMirroredBuf* metrics_buf_source,
                       GMEnv* env) {
  GMInsist(metrics && env);
  GMInsist(metrics_buf_target && metrics_buf_source);

  const int nvl = metrics->num_vector_local;

  const MPI_Datatype mpi_type = gm_mpi_type(env);

  const int mpi_code = MPI_Allreduce(metrics_buf_source->h,
                                     metrics_buf_target->h,
                                     nvl * (size_t)nvl, mpi_type, MPI_SUM,
                                     GMEnv_mpi_comm_field(env));
  GMInsist(mpi_code == MPI_SUCCESS);
}

//-----------------------------------------------------------------------------

MPI_Request gm_reduce_metrics_start(GMMetrics* metrics,
                                    GMMirroredBuf* metrics_buf_target,
                                    GMMirroredBuf* metrics_buf_source,
                                    GMEnv* env) {
  GMInsist(metrics && env);
  GMInsist(metrics_buf_target && metrics_buf_source);

  const int nvl = metrics->num_vector_local;

  const MPI_Datatype mpi_type = gm_mpi_type(env);

  MPI_Request mpi_request;
  const int mpi_code = MPI_Iallreduce(metrics_buf_source->h,
                                      metrics_buf_target->h,
                                      nvl * (size_t)nvl, mpi_type, MPI_SUM,
                                      GMEnv_mpi_comm_field(env),
                                      &mpi_request);
  GMInsist(mpi_code == MPI_SUCCESS);

  return mpi_request;
}

//-----------------------------------------------------------------------------

void gm_reduce_metrics_wait(MPI_Request* mpi_request, GMEnv* env) {
  GMInsist(mpi_request && env);

  MPI_Status mpi_status;

  const int mpi_code = MPI_Wait(mpi_request, &mpi_status);
  GMInsist(mpi_code == MPI_SUCCESS);
}

//=============================================================================
// Start/end transfer of vectors data to GPU

void gm_set_vectors_start(GMVectors* vectors, GMMirroredBuf* vectors_buf,
                          GMEnv* env) {
  GMInsist(vectors && vectors_buf && env);

  gm_linalg_set_matrix_start(vectors_buf, env);
}

//-----------------------------------------------------------------------------

void gm_set_vectors_wait(GMEnv* env) {
  GMInsist(env);

  gm_linalg_set_matrix_wait(env);
}

//=============================================================================
// Start/end transfer of metrics data from GPU

void gm_get_metrics_start(GMMetrics* metrics, GMMirroredBuf* metrics_buf,
                          GMEnv* env) {
  GMInsist(metrics && metrics_buf && env);

  gm_linalg_get_matrix_start(metrics_buf, env);
}

//-----------------------------------------------------------------------------

void gm_get_metrics_wait(GMMetrics* metrics, GMMirroredBuf* metrics_buf,
                         GMEnv* env) {
  GMInsist(metrics && metrics_buf && env);

  gm_linalg_get_matrix_wait(env);
}

//-----------------------------------------------------------------------------
