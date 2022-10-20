//-----------------------------------------------------------------------------
/*!
 * \file   decomp_mgr.hh
 * \author Wayne Joubert
 * \date   Tue Aug  8 19:58:57 EDT 2017
 * \brief  Define distribution of vectors to MPI ranks, padding needed, etc.
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

#include "cstdint"
#include "string.h"

#include "mpi.h"

#include "env.hh"
#include "tc.hh"
#include "decomp_mgr.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/// \brief Divisibility requirement for num_vector_local, due to sections.

size_t gm_nvl_divisibility_required_for_section(const CEnv& env) {

  const bool need_divisible_by_6 = env.num_way() == NumWay::_3 &&
                                   env.all2all() &&
                                   env.num_proc_vector() > 2;

  const size_t factor = need_divisible_by_6 ? 6 : 1;

  return factor;
}

//-----------------------------------------------------------------------------

size_t gm_nvl_size_required(size_t size_requested, const CEnv& env) {

  // NOTE: this function will in practice receive the same size_requested
  // on every rank and thus give the same result independent of MPI rank.

  // GEMMs on tensor cores sometimes impose divisibility conditions.
  const size_t factor_gemm = tc_nvl_divisibility_required_for_gemm(env);

  // For 3-way, the 6 sections of a block must all be the same size.
  const size_t factor_section = gm_nvl_divisibility_required_for_section(env);

  // Find LCM.  // (LATER: use std::lcm)

  COMET_ASSERT(factor_section == 1 || factor_section == 6);
  const size_t factor = factor_gemm % factor_section == 0 ? factor_gemm :
                        factor_gemm % 2 == 0 ? factor_gemm * 3 :
                        factor_gemm % 3 == 0 ? factor_gemm * 2 :
                                               factor_gemm * 6;

  const size_t size_for_gemm =
    tc_nvl_size_required_for_gemm(size_requested, env);

  // Pad up.
  return utils::ceil(size_for_gemm, factor) * factor;
}

//-----------------------------------------------------------------------------
// Set to null

GMDecompMgr GMDecompMgr_null() {
  GMDecompMgr result;
  memset((void*)&result, 0, sizeof(result));
  return result;
}

//-----------------------------------------------------------------------------
// (Pseudo) constructor

void GMDecompMgr_create(GMDecompMgr* dm,
                        bool fields_by_local,
                        bool vectors_by_local,
                        size_t num_field_specifier,
                        NV_t num_vector_specifier,
                        int vectors_data_type_id,
                        CEnv* env) {
  COMET_INSIST(dm && env);

  *dm = GMDecompMgr_null();

  if (!env->is_proc_active())
    return;

  //--------------------
  // Vector counts
  //--------------------

  if (vectors_by_local) {
    dm->num_vector_local = num_vector_specifier;
    const size_t nvl_required =
      gm_nvl_size_required(dm->num_vector_local, *env);
    COMET_INSIST_INTERFACE(env, dm->num_vector_local == nvl_required &&
         "Manual selection of nvl requires divisibility condition");
    // All vectors active on every proc.
    dm->num_vector_active_local = dm->num_vector_local;
    dm->num_vector_active = static_cast<NV_t>(dm->num_vector_active_local) *
                            env->num_proc_vector();
    dm->num_vector = static_cast<NV_t>(dm->num_vector_local) *
                     env->num_proc_vector();
    dm->vector_base = static_cast<NV_t>(dm->num_vector_local) *
                      env->proc_num_vector();
  } else { // ! vectors_by_local
    COMET_INSIST_INTERFACE(env, (env->all2all() ||
                                 env->num_proc_vector() == 1) &&
      "Use of all2all = no option currently requires "
      "num_vector_local, not num_vector.");
    dm->num_vector_active = num_vector_specifier;
    // Pad up as needed, require every proc has same number
    const int npv = env->num_proc_vector();
    const int pnv = env->proc_num_vector();
    //dm->num_vector_local = utils::ceil(dm->num_vector_active, (size_t)npv);
    dm->num_vector_local = gm_nvl_size_required(
      utils::ceil(dm->num_vector_active, static_cast<NV_t>(npv)), *env);
    dm->num_vector = static_cast<NV_t>(dm->num_vector_local) * npv;
    // Lower procs fully packed with active values
    // Upper procs fully inactive
    // One proc in between may be mixed
    const NV_t nvl = dm->num_vector_local;
    const NV_t nva = dm->num_vector_active;
    // Pack in lower procs with no gaps to ensure indexing of actives works
    // right independent of decomposition
    dm->num_vector_active_local = nva <= nvl * pnv ? 0 :
                                  nva >= nvl * (pnv + 1) ? nvl :
                                  nva - nvl * pnv;
    dm->vector_base = nvl * pnv <= nva ? nvl * pnv : nva;
    COMET_INSIST(nvl * pnv == dm->vector_base || 0 == dm->num_vector_active_local);
  } // if vectors_by_local

  // Notes on vector decompositon.
  //
  // If num_vector is specified (the typical production case) then:
  // nva = num_vector_active = the actual number requested by the user
  // nvl = num_vector_local = number of vectors stored per rank, to satisfy
  //   both nva and also divisibility conditions (same for all ranks).
  // nval = num_vector_active_local = number of active vectors on
  //   given rank (excluding padding for divisbility).
  //
  // The way it is set up now, the lower ranks are fully packed with active
  // vectors (nval = nvl) until vectors are exhausted, then possibly one
  // partly-full rank, then all ranks empty (all padding).
  // This is in some sense a simple arrangement.  However, it will be noted
  // different settings, e.g., tc, could require different padding values,
  // thus the parallel decomposition of vectors could be slightly different --
  // which would mean that partial runs based on stages and/or phases
  // could give different results.  Therefore, settings like "tc" and
  // "compute_method" should be kept constant across runs for such
  // run campaign, otherwise surprising results.
  //
  // This approach using padding is guaranteed to give exactly one copy
  // of each unique metric value, for the following reason:
  // 1) the original unpadded matrix (square) or tensor (cube) contains
  // a set of unique metrics, and additionally redundancies.
  // 2) this is mapped by a 1-1 map into a larger matrix or tensor,
  // by adding pad vectors at the high end.
  // 3) a subset of the entries of the larger (padded) tensor is taken to
  // eliminate redundancies.
  // 4) at the end, this is cropped back down to the original pre-padded
  // size.
  // 5) it is clear that no pair of elements of 4) can be redundant, since 
  // in such case so would 3).
  // 6) suppose some unique value of 4) is missing. Then also all of its
  // redundant copies must be missing.  This is to say, a tuple of 4)
  // is missing along with all permutations of its tuple indices.
  // But then they all must be missing from 3) since 3) only adds back
  // some tuples that have at least one index > nva. 3) must have all unique
  // (up to index ordering) tuples. Thus contradiction.

  //--------------------
  // Check the sizes
  //--------------------

  size_t sum = 0;

  COMET_INSIST(dm->num_vector_active >= 0 &&
           dm->num_vector_active <= dm->num_vector);
  COMET_INSIST(dm->num_vector_active_local >= 0 &&
           dm->num_vector_active_local <= dm->num_vector_local);

  COMET_MPI_SAFE_CALL(MPI_Allreduce(&dm->num_vector_local, &sum, 1,
    MPI_UNSIGNED_LONG_LONG, MPI_SUM, env->comm_repl_vector()));
  COMET_INSIST(sum == static_cast<NV_t>(dm->num_vector_local) *
                      env->num_proc_repl_vector() &&
           "Every process must have the same number of vectors.");
  COMET_INSIST(sum == dm->num_vector * env->num_proc_repl() &&
           "Error in local/global sizes computation.");

  COMET_MPI_SAFE_CALL(MPI_Allreduce(&dm->num_vector_active_local, &sum, 1,
    MPI_UNSIGNED_LONG_LONG, MPI_SUM, env->comm_repl_vector()));
  COMET_INSIST(sum == dm->num_vector_active * env->num_proc_repl() &&
           "Error in local/global sizes computation.");

  //--------------------
  // Field counts
  //--------------------

  if (fields_by_local) {
    dm->num_field_active_local = num_field_specifier;
    dm->num_field_local = dm->num_field_active_local;
    dm->num_field = dm->num_field_local * (size_t)env->num_proc_field();
    dm->num_field_active = dm->num_field;
    dm->field_base = dm->num_field_local * env->proc_num_field();
  } else { // ! fields_by_local
    dm->num_field_active = num_field_specifier;
    // Pad up as needed so that every proc has same number
    const int num_proc = env->num_proc_field();
    const int proc_num = env->proc_num_field();
    dm->num_field_local = utils::ceil(dm->num_field_active, (size_t)num_proc);
    dm->num_field = dm->num_field_local * num_proc;
    // Lower procs fully packed with active values
    // Upper procs fully inactive
    // One proc in between may be mixed
    // NOTE: see below for possible more inactive fields in packedfield
    const size_t nfl = dm->num_field_local;
    const size_t nfa = dm->num_field_active;
    dm->num_field_active_local =
      nfa <= nfl * proc_num       ? 0 :
      nfa >= nfl * (proc_num + 1) ? nfl :
                                    nfa - nfl * proc_num;
    // Pack in lower procs with no gaps to ensure indexing of actives works
    // right independent of decomposition
    dm->field_base = nfl * proc_num <= nfa ? nfl * proc_num : nfa;
    COMET_INSIST(nfl * proc_num == dm->field_base || 0 == dm->num_field_active_local);

  } // if fields_by_local

  //--------------------
  // Check the sizes
  //--------------------

  COMET_INSIST(dm->num_field_active >= 0 && dm->num_field_active <= dm->num_field);
  COMET_INSIST(dm->num_field_active_local >= 0 &&
           dm->num_field_active_local <= dm->num_field_local);

  COMET_MPI_SAFE_CALL(MPI_Allreduce(&dm->num_field_local, &sum, 1,
    MPI_UNSIGNED_LONG_LONG, MPI_SUM, env->comm_field()));
  COMET_INSIST(sum == dm->num_field_local * env->num_proc_field() &&
           "Every process must have the same number of fields.");
  COMET_INSIST(sum == dm->num_field &&
           "Error in local/global sizes computation.");

  COMET_MPI_SAFE_CALL(MPI_Allreduce(&dm->num_field_active_local, &sum, 1,
    MPI_UNSIGNED_LONG_LONG, MPI_SUM, env->comm_field()));
  COMET_INSIST(sum == dm->num_field_active &&
           "Error in local/global sizes computation.");

  //--------------------
  // Element sizes
  //--------------------

  switch (vectors_data_type_id) {
    //--------------------
    case DataTypeId::FLOAT: {
      dm->num_bit_per_field = BITS_PER_BYTE * sizeof(GMFloat);
      dm->num_bit_per_packedfield = BITS_PER_BYTE * sizeof(GMFloat);
    } break;
    //--------------------
    case DataTypeId::BITS2: {
      dm->num_bit_per_field = GM_BITS2_MAX_VALUE_BITS;
      dm->num_bit_per_packedfield = BITS_PER_BYTE * sizeof(GMBits2x64);
      // By design can only store this number of fields for this metric
      // TODO: later may be able to permit higher via rounding -
      // have 2-way method on-proc be exact, then for 3-way combining
      // or num_proc_field>1 drop low order bits to allow to fit.
      const int table_entry_value_max =
        env->metric_type() == MetricType::DUO ?  1 : env->pow2_num_way();
      COMET_INSIST_INTERFACE(env,
               ((uint64_t)(table_entry_value_max * dm->num_field)) <
                       (((uint64_t)1) << GM_TALLY1_MAX_VALUE_BITS)
                && "Number of fields requested is too large for this metric");
    } break;
    //--------------------
    default:
      COMET_INSIST(false && "Invalid vectors_data_type_id.");
  } //---switch---

  COMET_INSIST(dm->num_bit_per_packedfield % BITS_PER_BYTE == 0 &&
           "Error in size computation.");

  dm->num_field_per_packedfield = dm->num_bit_per_packedfield /
                                  dm->num_bit_per_field;

  //--------------------
  // Packedfield counts
  //--------------------

  const size_t d = tc_gemm_faxis_divisibility_required(*env);

  dm->num_packedfield_local =
    utils::ceil(utils::ceil(dm->num_field_local * dm->num_bit_per_field,
                            (size_t)dm->num_bit_per_packedfield), d) * d;

  //--------------------
  // Number of non-active fields on proc.
  //--------------------

  dm->num_pad_field_local =
    dm->num_packedfield_local *
    dm->num_field_per_packedfield -
    dm->num_field_active_local;

  //--------------------
  // tc memory
  //--------------------

  TCBufs::malloc(dm->num_vector_local,
                 dm->num_packedfield_local, dm->tc_bufs, *env);

  //--------------------
  // histograms
  //--------------------

  dm->histograms_default_ = new Histograms("", *env);

  dm->histograms_ = dm->histograms_default_;
}

//-----------------------------------------------------------------------------
// (Pseudo) destructor

void GMDecompMgr_destroy(GMDecompMgr* dm, CEnv* env) {
  COMET_INSIST(dm && env);

  if (! env->is_proc_active())
    return;

  //--------------------
  // tc memory
  //--------------------

  TCBufs::free(dm->tc_bufs, *env);

  //--------------------
  // histograms
  //--------------------

  dm->histograms_ = NULL;

  delete dm->histograms_default_;
  dm->histograms_default_ = NULL;
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
