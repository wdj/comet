//-----------------------------------------------------------------------------
/*!
 * \file   comm_xfer_utils.cc
 * \author Wayne Joubert
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Utilities for communication and transfer of vectors, metrics.
 */
//-----------------------------------------------------------------------------
/*-----------------------------------------------------------------------------

Copyright 2020, UT-Battelle, LLC

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

-----------------------------------------------------------------------------*/

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

void CommVectors::send_start(const Vectors& vectors,
                             int proc_num,
                             int mpi_tag,
                             CEnv& env) {
  COMET_INSIST(proc_num >= 0 && proc_num < env.num_proc_repl_vector());
  COMET_INSIST(!(env.is_comm_gpu() && !vectors.has_buf()));

  const int tag_multiplier = BuildHas::DEBUG ? 2 : 1;

  COMET_MPI_SAFE_CALL(MPI_Isend(
    env.is_comm_gpu() ?
      (void*)vectors.buf()->d :
      (void*)vectors.data(),
    vectors.num_packedfield_vector_local(), env.metrics_mpi_type(), proc_num,
    0+tag_multiplier*mpi_tag, env.comm_repl_vector(), &mpi_request_));

  if (BuildHas::DEBUG) {
    cksum_ = vectors.cksum();
    COMET_MPI_SAFE_CALL(MPI_Isend(
      (void*)&cksum_,
      sizeof(cksum_), MPI_BYTE, proc_num,
      1+2*mpi_tag, env.comm_repl_vector(), &mpi_request_cksum_));
  }
}

//-----------------------------------------------------------------------------
// Start/end MPI send/receive of vectors data

void CommVectors::recv_start(const Vectors& vectors,
                             int proc_num,
                             int mpi_tag,
                             CEnv& env) {
  COMET_INSIST(proc_num >= 0 && proc_num < env.num_proc_repl_vector());
  COMET_INSIST(!(env.is_comm_gpu() && !vectors.has_buf()));

  const int tag_multiplier = BuildHas::DEBUG ? 2 : 1;

  COMET_MPI_SAFE_CALL(MPI_Irecv(
    env.is_comm_gpu() ?
      (void*)vectors.buf()->d :
      (void*)vectors.data(),
    vectors.num_packedfield_vector_local(), env.metrics_mpi_type(), proc_num,
    0+tag_multiplier*mpi_tag, env.comm_repl_vector(), &mpi_request_));

  if (BuildHas::DEBUG) {
    COMET_MPI_SAFE_CALL(MPI_Irecv(
      (void*)&cksum_,
      sizeof(cksum_), MPI_BYTE, proc_num,
      1+2*mpi_tag, env.comm_repl_vector(), &mpi_request_cksum_));
  }
}

//=============================================================================
// MPI reduce operations

void reduce_metrics(GMMetrics& metrics,
                    MirroredBuf* target,
                    MirroredBuf* source,
                    CEnv& env) {
  COMET_INSIST(target && source);

  if (!env.do_reduce())
    return;

  const int nvl = metrics.num_vector_local;

  COMET_MPI_SAFE_CALL(MPI_Allreduce(source->h,
    target->h, nvl * (size_t)nvl, env.metrics_mpi_type(), MPI_SUM,
    env.comm_field()));
}

//-----------------------------------------------------------------------------

MPI_Request reduce_metrics_start(GMMetrics& metrics,
                                 MirroredBuf* target,
                                 MirroredBuf* source,
                                 CEnv& env) {
  COMET_INSIST(target && source);

  if (!env.do_reduce())
    return MPI_REQUEST_NULL;

  const int nvl = metrics.num_vector_local;

  target->lock_h();
  source->lock_h();

  MPI_Request mpi_request;
  COMET_MPI_SAFE_CALL(MPI_Iallreduce(source->h,
    target->h, nvl * (size_t)nvl, env.metrics_mpi_type(), MPI_SUM,
    env.comm_field(), &mpi_request));

  return mpi_request;
}

//-----------------------------------------------------------------------------

void reduce_metrics_wait(MPI_Request& mpi_request,
                         MirroredBuf* target,
                         MirroredBuf* source,
                         CEnv& env) {
  if (!env.do_reduce())
    return;

  MPI_Status mpi_status;

  COMET_MPI_SAFE_CALL(MPI_Wait(&mpi_request, &mpi_status));

  target->unlock_h();
  source->unlock_h();
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
