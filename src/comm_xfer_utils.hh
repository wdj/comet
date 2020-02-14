//-----------------------------------------------------------------------------
/*!
 * \file   comm_xfer_utils.hh
 * \author Wayne Joubert
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Utilities for communication and CPU/GPU transfer of vectors, metrics.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _comet_comm_xfer_utils_hh_
#define _comet_comm_xfer_utils_hh_

#include "env.hh"
#include "mirrored_buf.hh"
#include "vectors.hh"
#include "metrics.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------

MPI_Request gm_send_vectors_start(GMVectors* vectors, int proc_num, int mpi_tag,
                                  GMEnv* env);

MPI_Request gm_recv_vectors_start(GMVectors* vectors, int proc_num, int mpi_tag,
                                  GMEnv* env);

void gm_send_vectors_wait(MPI_Request* mpi_request, GMEnv* env);

void gm_recv_vectors_wait(MPI_Request* mpi_request, GMEnv* env);

//--------------------

void gm_reduce_metrics(GMMetrics* metrics,
                       MirroredBuf* metrics_buf_target,
                       MirroredBuf* metrics_buf_source,
                       GMEnv* env);

MPI_Request gm_reduce_metrics_start(GMMetrics* metrics,
                                    MirroredBuf* metrics_buf_target,
                                    MirroredBuf* metrics_buf_source,
                                    GMEnv* env);

void gm_reduce_metrics_wait(MPI_Request* mpi_request,
                            MirroredBuf* metrics_buf_target,
                            MirroredBuf* metrics_buf_source,
                            GMEnv* env);

//=============================================================================

} // namespace comet

#endif // _comet_comm_xfer_utils_hh_

//-----------------------------------------------------------------------------
