//-----------------------------------------------------------------------------
/*!
 * \file   decomp_mgr.hh
 * \author Wayne Joubert
 * \date   Tue Aug  8 19:58:57 EDT 2017
 * \brief  Manage parallel decomposition of vectors, fields, header.
 * \note   Copyright (C) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _gm_decomp_mgr_hh_
#define _gm_decomp_mgr_hh_

#include "env.hh"

//=============================================================================
// DecompMgr struct

typedef struct {
  // Field counts
  size_t num_field;
  size_t num_field_local;
  size_t num_field_active;
  size_t num_field_active_local;
  // Vector counts
  size_t num_vector;
  size_t num_vector_local;
  size_t num_vector_active;
  size_t num_vector_active_local;
  // Packed field info
  int num_bits_per_field;
  int num_bits_per_packedfield;
  int num_field_per_packedfield;
  int pad1;
  size_t num_packedfield_local;
} GMDecompMgr;

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

#endif // _gm_decomp_mgr_hh_

//-----------------------------------------------------------------------------
