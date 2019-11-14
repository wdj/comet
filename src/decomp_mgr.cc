//-----------------------------------------------------------------------------
/*!
 * \file   decomp_mgr.hh
 * \author Wayne Joubert
 * \date   Tue Aug  8 19:58:57 EDT 2017
 * \brief  Define distribution of vectors to MPI ranks, padding needed, etc.
 * \note   Copyright (C) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#include "cstdint"
#include "string.h"

#include "mpi.h"

#include "env.hh"
#include "linalg_tc.hh"
#include "decomp_mgr.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------

size_t gm_num_vector_local_required(size_t num_vector_active_local,
                                    GMEnv* const env) {
  GMInsist(env);
  // NOTE: this function should receive the same num_vector_active_local
  // and give the same result independent of MPI rank.

  const size_t factor_4 = gm_gemm_divisibility_required(env);

  const bool need_divisible_by_6 = env->num_way() == NUM_WAY::_3 &&
                                   env->all2all() &&
                                   env->num_proc_vector() > 2;

  const size_t lcm = (! need_divisible_by_6) ? factor_4 :
                     factor_4 % 2 == 0 ? 3 * factor_4 : 6 * factor_4;

  return gm_ceil_i8(num_vector_active_local, lcm)*lcm;
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
                        GMEnv* env) {
  GMInsist(dm && env);

  if (! env->is_proc_active()) {
    *dm = GMDecompMgr_null();
    return;
  }

  //const int vectors_data_type_id = env->data_type_vectors();

  //--------------------
  // Vector counts
  //--------------------

  if (vectors_by_local) {
    dm->num_vector_local = num_vector_specifier;
    const size_t num_vector_local_required = gm_num_vector_local_required(
                                              dm->num_vector_local, env);
    GMInsistInterface(env, dm->num_vector_local == num_vector_local_required &&
         "Manual selection of nvl requires divisibility condition");
    // All vectors active on every proc.
    dm->num_vector_active_local = dm->num_vector_local;
#if 0
    dm->num_vector_active_local = num_vector_specifier;
    dm->num_vector_local = gm_num_vector_local_required(
                               dm->num_vector_active_local, env);
#endif
    dm->num_vector_active = dm->num_vector_active_local *
                            env->num_proc_vector();
    dm->num_vector = dm->num_vector_local * env->num_proc_vector();
  } else { // ! vectors_by_local
    dm->num_vector_active = num_vector_specifier;
    // Pad up as needed, require every proc has same number
    const int num_proc = env->num_proc_vector();
    const int proc_num = env->proc_num_vector();
    //dm->num_vector_local = gm_ceil_i8(dm->num_vector_active, num_proc);
    dm->num_vector_local = gm_num_vector_local_required(
      gm_ceil_i8(dm->num_vector_active, env->num_proc_vector()), env);
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
  } // if vectors_by_local

  //--------------------
  // Check the sizes
  //--------------------

  size_t sum = 0;

  GMInsist(dm->num_vector_active >= 0 &&
           dm->num_vector_active <= dm->num_vector);
  GMInsist(dm->num_vector_active_local >= 0 &&
           dm->num_vector_active_local <= dm->num_vector_local);

  COMET_MPI_SAFE_CALL(MPI_Allreduce(&dm->num_vector_local, &sum, 1,
    MPI_UNSIGNED_LONG_LONG, MPI_SUM, env->comm_repl_vector()));
  GMInsist(sum == dm->num_vector_local * env->num_proc_repl_vector() &&
           "Every process must have the same number of vectors.");
  GMInsist(sum == dm->num_vector * env->num_proc_repl() &&
           "Error in local/global sizes computation.");

  COMET_MPI_SAFE_CALL(MPI_Allreduce(&dm->num_vector_active_local, &sum, 1,
    MPI_UNSIGNED_LONG_LONG, MPI_SUM, env->comm_repl_vector()));
  GMInsist(sum == dm->num_vector_active * env->num_proc_repl() &&
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
  } else { // ! fields_by_local
    dm->num_field_active = num_field_specifier;
    // Pad up as needed so that every proc has same number
    const int num_proc = env->num_proc_field();
    const int proc_num = env->proc_num_field();
    dm->num_field_local = gm_ceil_i8(dm->num_field_active, num_proc);
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

  } // if fields_by_local

  //--------------------
  // Check the sizes
  //--------------------

  GMInsist(dm->num_field_active >= 0 && dm->num_field_active <= dm->num_field);
  GMInsist(dm->num_field_active_local >= 0 &&
           dm->num_field_active_local <= dm->num_field_local);

  COMET_MPI_SAFE_CALL(MPI_Allreduce(&dm->num_field_local, &sum, 1,
    MPI_UNSIGNED_LONG_LONG, MPI_SUM, env->comm_field()));
  GMInsist(sum == dm->num_field_local * env->num_proc_field() &&
           "Every process must have the same number of fields.");
  GMInsist(sum == dm->num_field &&
           "Error in local/global sizes computation.");

  COMET_MPI_SAFE_CALL(MPI_Allreduce(&dm->num_field_active_local, &sum, 1,
    MPI_UNSIGNED_LONG_LONG, MPI_SUM, env->comm_field()));
  GMInsist(sum == dm->num_field_active &&
           "Error in local/global sizes computation.");

  //--------------------
  // Element sizes
  //--------------------

  const int bits_per_byte = 8;

  switch (vectors_data_type_id) {
    //--------------------
    case GM_DATA_TYPE_FLOAT: {
      dm->num_bits_per_field = bits_per_byte * sizeof(GMFloat);
      dm->num_bits_per_packedfield = bits_per_byte * sizeof(GMFloat);
    } break;
    //--------------------
    case GM_DATA_TYPE_BITS2: {
      dm->num_bits_per_field = GM_BITS2_MAX_VALUE_BITS;
      dm->num_bits_per_packedfield = bits_per_byte * sizeof(GMBits2x64);
      // By design can only store this number of fields for this metric
      // TODO: later may be able to permit higher via rounding -
      // have 2-way method on-proc be exact, then for 3-way combining
      // or num_proc_field>1 drop low order bits to allow to fit.
      const int table_entry_limit =
        env->metric_type() == MetricType::DUO ?
        1 : 1 << env->num_way();
      GMInsistInterface(env,
               ((uint64_t)(table_entry_limit * dm->num_field)) <
                       (((uint64_t)1) << GM_TALLY1_MAX_VALUE_BITS)
                && "Number of fields requested is too large for this metric");
    } break;
    //--------------------
    default:
      GMInsist(false && "Invalid vectors_data_type_id.");
  } //---switch---

  GMInsist(dm->num_bits_per_packedfield % bits_per_byte == 0 &&
           "Error in size computation.");

  dm->num_field_per_packedfield = dm->num_bits_per_packedfield /
                                  dm->num_bits_per_field;

  //--------------------
  // Packedfield counts
  //--------------------

  dm->num_packedfield_local =
      gm_ceil_i8(dm->num_field_local * dm->num_bits_per_field,
                 dm->num_bits_per_packedfield);

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

  gm_tc_bufs_malloc(dm->num_vector_local, dm->num_field_local,
                    dm->num_packedfield_local, dm->tc_bufs, env);
}

//-----------------------------------------------------------------------------
// (Pseudo) destructor

void GMDecompMgr_destroy(GMDecompMgr* dm, GMEnv* env) {
  GMInsist(dm && env);

  if (! env->is_proc_active()) {
    return;
  }

  //--------------------
  // tc memory
  //--------------------

  gm_tc_bufs_free(dm->tc_bufs, env);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
