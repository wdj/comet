//-----------------------------------------------------------------------------
/*!
 * \file   comm_xfer_utils.cc
 * \author Wayne Joubert
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Utilities for communication and CPU/GPU transfer of vectors, metrics.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#include "env.hh"
#include "linalg.hh"
#include "mirrored_buf.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "comm_xfer_utils.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
// Start/end MPI send/receive of vectors data

MPI_Request gm_send_vectors_start(GMVectors* vectors,
                                  int proc_num,
                                  int mpi_tag,
                                  GMEnv* env) {
  COMET_INSIST(vectors && env);
  COMET_INSIST(proc_num >= 0 && proc_num < env->num_proc_repl_vector());

  MPI_Request mpi_request;

  COMET_MPI_SAFE_CALL(MPI_Isend((void*)vectors->data,
    vectors->num_packedval_local, env->metrics_mpi_type(), proc_num, mpi_tag
  , env->comm_repl_vector(), &mpi_request));

  return mpi_request;
}

//-----------------------------------------------------------------------------

MPI_Request gm_recv_vectors_start(GMVectors* vectors,
                                  int proc_num,
                                  int mpi_tag,
                                  GMEnv* env) {
  COMET_INSIST(vectors && env);
  COMET_INSIST(proc_num >= 0 && proc_num < env->num_proc_repl_vector());

  MPI_Request mpi_request;

  COMET_MPI_SAFE_CALL(MPI_Irecv((void*)vectors->data,
    vectors->num_packedval_local, env->metrics_mpi_type(), proc_num, mpi_tag,
    env->comm_repl_vector(), &mpi_request));

  return mpi_request;
}

//-----------------------------------------------------------------------------

void gm_send_vectors_wait(MPI_Request* mpi_request, GMEnv* env) {
  COMET_INSIST(mpi_request && env);

  MPI_Status mpi_status;

  COMET_MPI_SAFE_CALL(MPI_Wait(mpi_request, &mpi_status));
}

//-----------------------------------------------------------------------------

void gm_recv_vectors_wait(MPI_Request* mpi_request, GMEnv* env) {
  COMET_INSIST(mpi_request && env);

  MPI_Status mpi_status;

  COMET_MPI_SAFE_CALL(MPI_Wait(mpi_request, &mpi_status));
}

//=============================================================================
// MPI reduce operations

void gm_reduce_metrics(GMMetrics* metrics,
                       MirroredBuf* metrics_buf_target,
                       MirroredBuf* metrics_buf_source,
                       GMEnv* env) {
  COMET_INSIST(metrics && env);
  COMET_INSIST(metrics_buf_target && metrics_buf_source);

  if (!env->do_reduce())
    return;

  const int nvl = metrics->num_vector_local;

  COMET_MPI_SAFE_CALL(MPI_Allreduce(metrics_buf_source->h,
    metrics_buf_target->h, nvl * (size_t)nvl, env->metrics_mpi_type(), MPI_SUM,
    env->comm_field()));
}

//-----------------------------------------------------------------------------

MPI_Request gm_reduce_metrics_start(GMMetrics* metrics,
                                    MirroredBuf* metrics_buf_target,
                                    MirroredBuf* metrics_buf_source,
                                    GMEnv* env) {
  COMET_INSIST(metrics && env);
  COMET_INSIST(metrics_buf_target && metrics_buf_source);

  if (!env->do_reduce())
    return MPI_REQUEST_NULL;

  const int nvl = metrics->num_vector_local;

  metrics_buf_target->lock_h();
  metrics_buf_source->lock_h();

  MPI_Request mpi_request;
  COMET_MPI_SAFE_CALL(MPI_Iallreduce(metrics_buf_source->h,
    metrics_buf_target->h, nvl * (size_t)nvl, env->metrics_mpi_type(), MPI_SUM,
    env->comm_field(), &mpi_request));

  return mpi_request;
}

//-----------------------------------------------------------------------------

void gm_reduce_metrics_wait(MPI_Request* mpi_request,
                            MirroredBuf* metrics_buf_target,
                            MirroredBuf* metrics_buf_source,
                            GMEnv* env) {
  COMET_INSIST(mpi_request && env);

  if (!env->do_reduce())
    return;

  MPI_Status mpi_status;

  COMET_MPI_SAFE_CALL(MPI_Wait(mpi_request, &mpi_status));

  metrics_buf_target->unlock_h();
  metrics_buf_source->unlock_h();
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
