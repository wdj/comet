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

  if(env.print_details()) printf("factor_gemm=%zu factor_section=%zu factor=%zu size_for_gemm=%zu size_req=%zu\n",
    factor_gemm,factor_section,factor,size_for_gemm,
    utils::ceil(size_for_gemm, factor)*factor);
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
                        size_t num_vector_specifier,
                        int vectors_data_type_id,
                        CEnv* env) {
  COMET_INSIST(dm && env);

  if (! env->is_proc_active()) {
    *dm = GMDecompMgr_null();
    return;
  }

  if(env->print_details()) printf("In DecompManager vectors_by_local=%d fields_by_local=%d\n",vectors_by_local,fields_by_local);

  //--------------------
  // Vector counts
  //--------------------

  if (vectors_by_local) {
    dm->num_vector_local = num_vector_specifier;
    const size_t nvl_required =
      gm_nvl_size_required(dm->num_vector_local, *env);
    if(env->print_details()) printf("In vectors_by_local num_vector_local=%zu nvl_required=%zu\n",
      dm->num_vector_local,nvl_required);
    COMET_INSIST_INTERFACE(env, dm->num_vector_local == nvl_required &&
         "Manual selection of nvl requires divisibility condition");
    // All vectors active on every proc.
    dm->num_vector_active_local = dm->num_vector_local;
    dm->num_vector_active = dm->num_vector_active_local *
                            env->num_proc_vector();
    dm->num_vector = dm->num_vector_local * env->num_proc_vector();
    dm->vector_base = dm->num_vector_local * env->proc_num_vector();
    if(env->print_details())
      printf("In vectors_by_local num_vector=%zu num_vector_active=%zu\n",
             dm->num_vector,dm->num_vector_active);
  } else { // ! vectors_by_local
    COMET_INSIST_INTERFACE(env, (env->all2all() ||
                                 env->num_proc_vector() == 1) &&
      "Use of all2all = no option currently requires "
      "num_vector_local, not num_vector.");
    dm->num_vector_active = num_vector_specifier;
    // Pad up as needed, require every proc has same number
    const int num_proc = env->num_proc_vector();
    const int proc_num = env->proc_num_vector();
    //dm->num_vector_local = utils::ceil(dm->num_vector_active, (size_t)num_proc);
    dm->num_vector_local = gm_nvl_size_required(
      utils::ceil(dm->num_vector_active, (size_t)env->num_proc_vector()), *env);
    dm->num_vector = dm->num_vector_local * num_proc;
    // Lower procs fully packed with active values
    // Upper procs fully inactive
    // One proc in between may be mixed
    const size_t nvl = dm->num_vector_local;
    const size_t nva = dm->num_vector_active;
    // Pack in lower procs with no gaps to ensure indexing of actives works
    // right independent of decomposition
    dm->num_vector_active_local = nva <= nvl * proc_num ? 0 :
                                  nva >= nvl * (proc_num + 1) ? nvl :
                                  nva - nvl * proc_num;
    dm->vector_base = nvl * proc_num <= nva ? nvl * proc_num : nva;
    COMET_INSIST(nvl * proc_num == dm->vector_base || 0 == dm->num_vector_active_local);
    if(env->print_details()) printf("In !vectors_by_local num_vector_active=%zu num_vector_local=%zu num_vector=%zu num_vector_active_local=%zu vector_base=%zu\n",dm->num_vector_active,dm->num_vector_local,dm->num_vector,dm->num_vector_active_local,dm->vector_base);
  } // if vectors_by_local

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
  COMET_INSIST(sum == dm->num_vector_local * env->num_proc_repl_vector() &&
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
    dm->num_field = dm->num_field_local * env->num_proc_field();
    dm->num_field_active = dm->num_field;
    dm->field_base = dm->num_field_local * env->proc_num_field();
    if(env->print_details())
      printf("In fields_by_local num_field=%zu num_field_local=%zu\n",
             dm->num_field,dm->num_field_local);
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
    if(env->print_details()) printf("In !fields_by_local\n");
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

  const int bits_per_byte = 8;

  switch (vectors_data_type_id) {
    //--------------------
    case GM_DATA_TYPE_FLOAT: {
      dm->num_bit_per_field = bits_per_byte * sizeof(GMFloat);
      dm->num_bit_per_packedfield = bits_per_byte * sizeof(GMFloat);
      if(env->print_details())
        printf("In GM_DATA_TYPE_FLOAT bits_per_byte=%d size=%zu\n",
               bits_per_byte,sizeof(GMFloat));
    } break;
    //--------------------
    case GM_DATA_TYPE_BITS2: {
      dm->num_bit_per_field = GM_BITS2_MAX_VALUE_BITS;
      dm->num_bit_per_packedfield = bits_per_byte * sizeof(GMBits2x64);
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
      if(env->print_details())
        printf("In GM_DATA_TYPE_BITS2 num_bit_per_field=%d bits_per_byte=%d "
               "GMBits2x64=%zu num_bit_per_packedfield=%d\n",
               dm->num_bit_per_field,bits_per_byte,sizeof(GMBits2x64),dm->num_bit_per_packedfield);
    } break;
    //--------------------
    default:
      COMET_INSIST(false && "Invalid vectors_data_type_id.");
  } //---switch---

  COMET_INSIST(dm->num_bit_per_packedfield % bits_per_byte == 0 &&
           "Error in size computation.");

  dm->num_field_per_packedfield = dm->num_bit_per_packedfield /
                                  dm->num_bit_per_field;
  if(env->print_details())
    printf("num_field_per_packedfield = %d/%d = %d\n",
           dm->num_bit_per_packedfield,dm->num_bit_per_field,dm->num_field_per_packedfield);

  //--------------------
  // Packedfield counts
  //--------------------

  const size_t d = tc_gemm_faxis_divisibility_required(*env);

  dm->num_packedfield_local =
    utils::ceil(utils::ceil(dm->num_field_local * dm->num_bit_per_field,
                            (size_t)dm->num_bit_per_packedfield), d) * d;
  if(env->print_details())
    printf("num_packedfield_local = ceil(%zu*%d,%d)=%zu\n",dm->num_field_local,
           dm->num_bit_per_field, dm->num_bit_per_packedfield,dm->num_packedfield_local);

  //--------------------
  // Number of non-active fields on proc.
  //--------------------

  dm->num_pad_field_local =
    dm->num_packedfield_local *
    dm->num_field_per_packedfield -
    dm->num_field_active_local;
  if(env->print_details()) 
    printf("num_pad_field_local = %zu*%d - %zu = %d\n",dm->num_packedfield_local,
           dm->num_field_per_packedfield,dm->num_field_active_local,dm->num_pad_field_local);

  //--------------------
  // tc memory
  //--------------------
  if(env->print_details()) printf("Calling TCBufs malloc\n");
  TCBufs::malloc(dm->num_vector_local, //dm->num_field_local,
                 dm->num_packedfield_local, dm->tc_bufs, *env);

  if(env->print_details()) printf("Done in DecompManager\n");
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
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
