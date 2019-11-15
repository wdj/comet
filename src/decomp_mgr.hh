//-----------------------------------------------------------------------------
/*!
 * \file   decomp_mgr.hh
 * \author Wayne Joubert
 * \date   Tue Aug  8 19:58:57 EDT 2017
 * \brief  Define distribution of vectors to MPI ranks, padding needed, etc.
 * \note   Copyright (C) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _comet_decomp_mgr_hh_
#define _comet_decomp_mgr_hh_

#include "env.hh"
#include "linalg_tc.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
// DecompMgr struct

typedef struct {
  // Field counts
  size_t num_field;
  size_t num_field_local;
  size_t num_field_active;
  size_t num_field_active_local;
  size_t field_base;
  // Vector counts
  size_t num_vector;
  size_t num_vector_local;
  size_t num_vector_active;
  size_t num_vector_active_local;
  // Packed field info
  int num_bits_per_field;
  int num_bits_per_packedfield;
  int num_field_per_packedfield;
  int num_pad_field_local;
  size_t num_packedfield_local;
  // Other
  TCBufs tc_bufs;
} GMDecompMgr;

//-----------------------------------------------------------------------------

size_t gm_num_vector_local_required(size_t num_vector_active,
                                    GMEnv* const env);

//-----------------------------------------------------------------------------
// Set to null

GMDecompMgr GMDecompMgr_null();

//-----------------------------------------------------------------------------
// (Pseudo) constructor

void GMDecompMgr_create(GMDecompMgr* dm,
                        bool fields_by_local,
                        bool vectors_by_local,
                        size_t num_field_specifier,
                        size_t num_vector_specifier,
                        int vectors_data_type_id,
                        GMEnv* env);

//-----------------------------------------------------------------------------
// (Pseudo) destructor

void GMDecompMgr_destroy(GMDecompMgr* dm, GMEnv* env);

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _comet_decomp_mgr_hh_

//-----------------------------------------------------------------------------
