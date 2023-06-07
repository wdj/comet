//-----------------------------------------------------------------------------
/*!
 * \file   comm_xfer_utils.hh
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

#ifndef _COMET_COMM_XFER_UTILS_HH_
#define _COMET_COMM_XFER_UTILS_HH_

#include "env.hh"
#include "mirrored_buf.hh"
#include "vectors.hh"
#include "metrics.hh"

//-----------------------------------------------------------------------------

namespace comet {

//=============================================================================

class CommVectors {

public:

  CommVectors()
    : mpi_request_(MPI_REQUEST_NULL)
    , mpi_request_cksum_(MPI_REQUEST_NULL)
    , cksum_(0) {
  }

  void wait() {
    MPI_Status mpi_status;
    COMET_MPI_SAFE_CALL(MPI_Wait(&mpi_request_, &mpi_status));
    if (BuildHas::DEBUG) {
      MPI_Status mpi_status_cksum;
      COMET_MPI_SAFE_CALL(MPI_Wait(&mpi_request_cksum_, &mpi_status_cksum));
    }
  }

  void send_start(const Vectors& vectors,
                  int proc_num,
                  int mpi_tag,
                  CEnv& env);

  void recv_start(const Vectors& vectors,
                  int proc_num,
                  int mpi_tag,
                  CEnv& env);

  size_t cksum() const {return cksum_;}

private:

  MPI_Request mpi_request_;
  MPI_Request mpi_request_cksum_;
  size_t cksum_;

}; // CommVectors

//-----------------------------------------------------------------------------

void reduce_metrics(MirroredBuf* target,
                    MirroredBuf* source,
                    CEnv& env);

MPI_Request reduce_metrics_start(MirroredBuf* target,
                                 MirroredBuf* source,
                                 CEnv& env);

void reduce_metrics_wait(MPI_Request& mpi_request,
                         MirroredBuf* target,
                         MirroredBuf* source,
                         CEnv& env);

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_COMM_XFER_UTILS_HH_

//-----------------------------------------------------------------------------
