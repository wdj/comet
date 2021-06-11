//-----------------------------------------------------------------------------
/*!
 * \file   histograms.cc
 * \author Wayne Joubert
 * \date   Fri Jun 11 08:15:19 EDT 2021
 * \brief  Manage histograms, definitions
 */
//-----------------------------------------------------------------------------
/*-----------------------------------------------------------------------------

Copyright 2021, UT-Battelle, LLC

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
#include "mirrored_buf.hh"
#include "histograms.hh"

//-----------------------------------------------------------------------------

namespace comet {

//-----------------------------------------------------------------------------
/// \brief Constructor for Histograms class.

Histograms::Histograms(CEnv& env)
  : env_(env)
  , range_(env_.ccc_duo_multiplier())
  , num_buckets_((int)(RECIP_BUCKET_WIDTH * range_))
  , buf_(num_buckets_, num_histograms(), env_) {

  COMET_INSIST(env_.is_metric_type_bitwise());




}

//-----------------------------------------------------------------------------

void Histograms::output() {

}

//-----------------------------------------------------------------------------

void Histograms::reduce() {

}



} // namespace comet

//=============================================================================

//-----------------------------------------------------------------------------

