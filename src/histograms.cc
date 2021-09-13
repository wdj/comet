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

Histograms::Histograms(const char* histograms_file, CEnv& env)
  : env_(env)
  , is_computing_histograms_(histograms_file && env_.is_metric_type_bitwise())
  , histograms_file_str_(histograms_file ? histograms_file : "")
  , range_(is_computing_histograms_ ? env_.ccc_duo_multiplier() : 0)
  , num_buckets_((int)(RECIP_BUCKET_WIDTH * range_))
  , num_histograms_(env_.num_way() + 2)
  , num_elts_(num_buckets_ * num_histograms_)
  , buf_(num_buckets_, num_histograms(), sizeof(Elt_t),  env_)
  , buf_finalized_(num_buckets_, num_histograms(), sizeof(Elt_t),  env_)
  , is_finalized_(false) {

  //COMET_INSIST(env_.is_metric_type_bitwise() || !histograms_file);

  if (!is_computing_histograms() || !env_.is_proc_active())
    return;

  buf_.set_zero_h();
  buf_.to_accel();
}

//-----------------------------------------------------------------------------
/// \brief Finalize computation of histograms.

void Histograms::finalize() {

  if (!is_computing_histograms() || !env_.is_proc_active())
    return;

  if (is_finalized_)
    return;

  // Retrieve from accelerator if necessary.

  if (is_computing_on_accel())
    buf_.from_accel();

  // MPI reduce of histograms across proc_repl_vector ranks.

  COMET_INSIST(sizeof(Elt_t) == sizeof(double));

  COMET_MPI_SAFE_CALL(MPI_Allreduce(buf_.h, buf_finalized_.h,
    num_elts_, MPI_DOUBLE, MPI_SUM, env_.comm_repl_vector()));

  is_finalized_ = true;
}

//-----------------------------------------------------------------------------
/// \brief Write histograms to file.

void Histograms::output() {

  if (!is_computing_histograms() || !env_.is_proc_active())
    return;

  COMET_INSIST(is_finalized_);

  if (env_.proc_num() != 0)
    return;

  FILE* file = fopen(histograms_file_str_.c_str(), "w");
  COMET_INSIST(NULL != file && "Unable to open file.");

  // Write file header.
  if (env_.num_way() == 2)
    fprintf(file, "min	LL	LH	HH	LL+HH\n");
  else
    fprintf(file, "min	LLL	LLH	LHH	HHH	LLL+HHH\n");

  // Write line for each hist. bucket.
  for (int row = 0; row < num_buckets_ ; ++row) {
    fprintf(file, "%.3f	",  bucket_min(row));
    for (int col = 0; col < num_histograms(); ++col) {
      fprintf(file, "%.0f%s", elt_finalized(row, col),
        col < num_histograms()-1 ? "\t" : "\n");
    }
  }

  fclose(file);
}

//-----------------------------------------------------------------------------
/// \brief Check result.

void Histograms::check(size_t num_vector_active) {

  if (!is_computing_histograms() || !env_.is_proc_active())
    return;

  COMET_INSIST(is_finalized_);

  if (env_.proc_num() != 0)
    return;

  // NOTE: this will only work if computing all stages and phases.

  const size_t num_metrics_total = env_.num_way() == 2 ?
    ((num_vector_active) * (num_vector_active - 1))  / 2 :
    ((num_vector_active) * (num_vector_active - 1) *
     (num_vector_active - 2)) / 6;

  const size_t expected = num_metrics_total * (1 << env_.num_way());
  const size_t calculated = sum_();

  if (calculated != expected)
    fprintf(stderr, "Error in histogram calculation: "
      "expected %zu, calculated %zu.\n", expected, calculated);

  COMET_INSIST(calculated == expected);
}

//-----------------------------------------------------------------------------
/// \brief Compute sum of all (finalized) histogram elements.

Histograms::Elt_t Histograms::sum_() const {
  COMET_INSIST(is_finalized_);
  Elt_t result = 0;
  for (int col = 0; col < num_histograms_; ++col) {
    const int multiplier = env_.num_way() + 1 == col ? 0 : 1;
    for (int row = 0; row < num_buckets_ ; ++row) {
      result += multiplier * elt_finalized_const(row, col);
    }
  }
  return result;
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

