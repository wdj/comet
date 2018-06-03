//-----------------------------------------------------------------------------
/*!
 * \file   decomp_mgr.hh
 * \author Wayne Joubert
 * \date   Tue Aug  8 19:58:57 EDT 2017
 * \brief  Manage parallel decomposition of vectors, fields, header.
 * \note   Copyright (C) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#include "string.h"

#include "mpi.h"

#include "env.hh"
#include "decomp_mgr.hh"

//=============================================================================
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

  if (! GMEnv_is_proc_active(env)) {
    *dm = GMDecompMgr_null();
    return;
  }

  //const int vectors_data_type_id = GMEnv_data_type_vectors(env);

  //--------------------
  // Vector counts
  //--------------------

  if (vectors_by_local) {
    dm->num_vector_active_local = num_vector_specifier;
    dm->num_vector_local = dm->num_vector_active_local;
    dm->num_vector = dm->num_vector_local * GMEnv_num_proc_vector_i(env);
    dm->num_vector_active = dm->num_vector;
  } else { // ! vectors_by_local
    dm->num_vector_active = num_vector_specifier;
    // Pad up as needed, require every proc has same number
    const int num_proc = GMEnv_num_proc_vector_i(env);
    const int proc_num = GMEnv_proc_num_vector_i(env);
    //dm->num_vector_local = gm_ceil_i8(dm->num_vector_active, num_proc);
    dm->num_vector_local = gm_num_vector_local_required(dm->num_vector_active,
                                                        env);
    dm->num_vector = dm->num_vector_local * num_proc;
    // Lower procs fully packed with active values
    // Upper procs fully inactive
    // One proc in between may be mixed
    const size_t nvl = dm->num_vector_local;
    const size_t nva = dm->num_vector_active;
    dm->num_vector_active_local = nva <= nvl * proc_num ? 0 :
                                  nva >= nvl * (proc_num + 1) ? nvl :
                                  nva - nvl * proc_num;
  } // if vectors_by_local

  //--------------------
  // Check the sizes
  //--------------------

  int mpi_code = 0;
  size_t sum = 0;

  GMInsist(dm->num_vector_active >= 0 &&
           dm->num_vector_active <= dm->num_vector);
  GMInsist(dm->num_vector_active_local >= 0 &&
           dm->num_vector_active_local <= dm->num_vector_local);

  mpi_code = MPI_Allreduce(&dm->num_vector_local, &sum, 1,
                           MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                           GMEnv_mpi_comm_repl_vector(env));
  GMInsist(mpi_code == MPI_SUCCESS);
  GMInsist(sum == dm->num_vector_local * GMEnv_num_proc_repl_vector(env) &&
           "Every process must have the same number of vectors.");
  GMInsist(sum == dm->num_vector * GMEnv_num_proc_repl(env));

  mpi_code = MPI_Allreduce(&dm->num_vector_active_local, &sum, 1,
                           MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                           GMEnv_mpi_comm_repl_vector(env));
  GMInsist(mpi_code == MPI_SUCCESS);
  GMInsist(sum == dm->num_vector_active * GMEnv_num_proc_repl(env));

  //--------------------
  // Field counts
  //--------------------

  if (fields_by_local) {
    dm->num_field_active_local = num_field_specifier;
    dm->num_field_local = dm->num_field_active_local;
    dm->num_field = dm->num_field_local * GMEnv_num_proc_field(env);
    dm->num_field_active = dm->num_field;
    dm->field_base = dm->num_field_local * GMEnv_proc_num_field(env);
  } else { // ! fields_by_local
    dm->num_field_active = num_field_specifier;
    // Pad up as needed so that every proc has same number
    const int num_proc = GMEnv_num_proc_field(env);
    const int proc_num = GMEnv_proc_num_field(env);
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
    dm->field_base = nfl * proc_num <= nfa ? nfl * proc_num : nfa;

  } // if fields_by_local

  //--------------------
  // Check the sizes
  //--------------------

  GMInsist(dm->num_field_active >= 0 && dm->num_field_active <= dm->num_field);
  GMInsist(dm->num_field_active_local >= 0 &&
           dm->num_field_active_local <= dm->num_field_local);

  mpi_code = MPI_Allreduce(&dm->num_field_local, &sum, 1,
                           MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                           GMEnv_mpi_comm_field(env));
  GMInsist(mpi_code == MPI_SUCCESS);
  GMInsist(sum == dm->num_field_local * GMEnv_num_proc_field(env) &&
           "Every process must have the same number of fields.");
  GMInsist(sum == dm->num_field);

  mpi_code = MPI_Allreduce(&dm->num_field_active_local, &sum, 1,
                           MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                           GMEnv_mpi_comm_field(env));
  GMInsist(mpi_code == MPI_SUCCESS);
  GMInsist(sum == dm->num_field_active);

  //--------------------
  // Element sizes
  //--------------------

  const int bits_per_byte = 8;

  switch (vectors_data_type_id) {
    case GM_DATA_TYPE_FLOAT: {
      dm->num_bits_per_field = bits_per_byte * sizeof(GMFloat);
      dm->num_bits_per_packedfield = bits_per_byte * sizeof(GMFloat);
    } break;
    case GM_DATA_TYPE_BITS2: {
      dm->num_bits_per_field = GM_BITS2_MAX_VALUE_BITS;
      dm->num_bits_per_packedfield = bits_per_byte * sizeof(GMBits2x64);
      /*---By design can only store this number of fields for this metric---*/
      GMInsistInterface(env,
               ((GMUInt64)(4 * dm->num_field)) <
                       (((GMUInt64)1) << GM_TALLY1_MAX_VALUE_BITS)
                && "Number of fields requested is too large for this metric");
    } break;
    /*--------------------*/
    default:
      GMInsist(false && "Invalid data type.");
  } /*---switch---*/

  GMInsist(dm->num_bits_per_packedfield % bits_per_byte == 0);

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
}

//-----------------------------------------------------------------------------
// (Pseudo) destructor

void GMDecompMgr_destroy(GMDecompMgr* dm, GMEnv* env) {
  GMInsist(dm && env);

  if (! GMEnv_is_proc_active(env)) {
    return;
  }

}

//-----------------------------------------------------------------------------
