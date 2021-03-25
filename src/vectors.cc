//-----------------------------------------------------------------------------
/*!
 * \file   vectors.cc
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Class to manage the set of vectors taken as input to the methods.
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

#include "cstdlib"
#include "cstdint"
#include "cstdio"
#include "string.h"

#include "mpi.h"

#include "env.hh"
#include "decomp_mgr.hh"
#include "vectors.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
// Null object.

GMVectors GMVectors_null() {
  GMVectors result;
  memset((void*)&result, 0, sizeof(GMVectors));
  return result;
}

//=============================================================================
// Set vector entries to zero.

void GMVectors_initialize(GMVectors* vectors, CEnv* env) {
  COMET_INSIST(vectors && env);

  if (! env->is_proc_active()) {
    return;
  }

  const size_t pfl_min = 0;
  const size_t pfl_max = vectors->dm->num_packedfield_local;

  switch (vectors->data_type_id) {
    case GM_DATA_TYPE_FLOAT: {
      GMFloat zero = 0;
      for (int vl = 0; vl < vectors->num_vector_local; ++vl) {
        for (size_t pfl = pfl_min; pfl < pfl_max; ++pfl) {
          GMVectors_float_set(vectors, pfl, vl, zero, env);
        }
      }
    } break;
    case GM_DATA_TYPE_BITS2: {
      const GMBits2x64 zero = GMBits2x64_null();
      for (int vl = 0; vl < vectors->num_vector_local; ++vl) {
        for (size_t pfl = pfl_min; pfl < pfl_max; ++pfl) {
          GMVectors_bits2x64_set(vectors, pfl, vl, zero, env);
        } // for pfl
      }
    } break;
    default:
      COMET_INSIST(false && "Invalid vectors data_type_id.");
  } // switch
}

//=============================================================================
// Set unused (pad) vector entries to zero.

//TODO: initialize pad vectors as well as pad fields !!!

// TODO: consider check to make sure user has set all vector entries (?).

void GMVectors_initialize_pad(GMVectors* vectors, CEnv* env) {
  COMET_INSIST(vectors && env);

  /*---Ensure final pad words/bits of each vector are set to zero so that
       word-wise summations of bits aren't corrupted with bad trailing data---*/

  if (! env->is_proc_active()) {
    return;
  }

  const size_t nfal = vectors->dm->num_field_active_local;

  const size_t pfl_min = nfal / vectors->dm->num_field_per_packedfield;
  const size_t pfl_max = vectors->dm->num_packedfield_local;

  switch (vectors->data_type_id) {
    case GM_DATA_TYPE_FLOAT: {
      GMFloat zero = 0;
      for (int vl = 0; vl < vectors->num_vector_local; ++vl) {
        for (size_t pfl = pfl_min; pfl < pfl_max; ++pfl) {
          GMVectors_float_set(vectors, pfl, vl, zero, env);
        }
      }
    } break;
    case GM_DATA_TYPE_BITS2: {
      const GMBits2x64 zero = GMBits2x64_null();
      const uint64_t allbits = 0xffffffffffffffff;
      for (int vl = 0; vl < vectors->num_vector_local; ++vl) {
        for (size_t pfl = pfl_min; pfl < pfl_max; ++pfl) {
          const size_t fl = 64 * pfl;

          if (fl >= nfal) {

            GMVectors_bits2x64_set(vectors, pfl, vl, zero, env);

          } else if (fl + 32 >= nfal) {

            GMBits2x64 val = GMVectors_bits2x64_get(vectors, pfl, vl, env);
            const int shift_dist = 64 - 2*(nfal-fl);
            COMET_ASSERT(shift_dist >= 0 && shift_dist < 64);
            val.data[0] &= allbits >> shift_dist;
            val.data[1] = 0;
            GMVectors_bits2x64_set(vectors, pfl, vl, val, env);

          } else if (fl + 64 >= nfal) {

            GMBits2x64 val = GMVectors_bits2x64_get(vectors, pfl, vl, env);
            const int shift_dist = 64 - 2*(nfal-fl-32);
            COMET_ASSERT(shift_dist >= 0 && shift_dist < 64);
            val.data[1] &= allbits >> shift_dist;
            GMVectors_bits2x64_set(vectors, pfl, vl, val, env);

          } // if
        } // for pfl
      }
    } break;
    default:
      COMET_INSIST(false && "Invalid vectors data_type_id.");
  } // switch
}

//=============================================================================
// Vectors pseudo-constructor.

void GMVectors_create_imp_(GMVectors* vectors,
                           int data_type_id,
                           GMDecompMgr* dm,
                           CEnv* env) {
  COMET_INSIST(vectors && dm && env);

  vectors->data_type_id = data_type_id;
  vectors->dm = dm;
  vectors->num_field = dm->num_field;
  vectors->num_field_active = dm->num_field_active;
  vectors->num_field_local = dm->num_field_local;
  vectors->num_vector = dm->num_vector;
  vectors->num_vector_local = dm->num_vector_local;

  const int bits_per_byte = 8;

  // Allocation size for vector storage

  vectors->num_packedfield_local = dm->num_packedfield_local;

  vectors->num_packedfield_vector_local =
      vectors->num_packedfield_local * dm->num_vector_local;

  vectors->data_size = vectors->num_packedfield_vector_local *
                       (dm->num_bit_per_packedfield / bits_per_byte);

  // Set up vector storage, mirrored buffer

  vectors->buf_ = new MirroredBuf(*env);

  if (vectors->has_buf_) {
    vectors->buf_->allocate(vectors->num_packedfield_local,
                            vectors->num_vector_local);
    vectors->data = vectors->buf_->h; // alias vector storage to buf
  } else {
    vectors->data = gm_malloc(vectors->data_size, env);
  }

  // Set pad entries to zero

  GMVectors_initialize_pad(vectors, env);
}

//-----------------------------------------------------------------------------

void GMVectors_create(GMVectors* vectors,
                      int data_type_id,
                      GMDecompMgr* dm,
                      CEnv* env) {
  COMET_INSIST(vectors && dm && env);

  *vectors = GMVectors_null();

  if (! env->is_proc_active()) {
    return;
  }

  vectors->has_buf_ = false;

  GMVectors_create_imp_(vectors, data_type_id, dm, env);
}

//-----------------------------------------------------------------------------

void GMVectors_create_with_buf(GMVectors* vectors,
                               int data_type_id,
                               GMDecompMgr* dm,
                               CEnv* env) {
  COMET_INSIST(vectors && dm && env);

  *vectors = GMVectors_null();

  if (! env->is_proc_active()) {
    return;
  }

  //vectors->has_buf_ = env->is_using_linalg();
  vectors->has_buf_ = true;

  GMVectors_create_imp_(vectors, data_type_id, dm, env);
}

//=============================================================================
// Vectors pseudo-destructor.

void GMVectors_destroy(GMVectors* vectors, CEnv* env) {
  COMET_INSIST(vectors && env);
  COMET_INSIST(vectors->data || ! env->is_proc_active());

  if (! env->is_proc_active()) {
    return;
  }

  if (!vectors->has_buf_) {
    gm_free(vectors->data, vectors->data_size, env);
  }

  delete vectors->buf_;

  *vectors = GMVectors_null();
}


//=============================================================================
// Copy vectors to mirrored buffer

void gm_vectors_to_buf(MirroredBuf* vectors_buf,
                       GMVectors* vectors,
                       CEnv* env) {
  COMET_INSIST(vectors && vectors_buf && env);

  if (!env->is_using_linalg())
    return;

  // Copy vectors into GPU buffers if needed.

  switch (env->metric_type()) {
    case MetricType::CZEK: {
      // don't use collapse because of overflow for large sizes
      //#pragma omp parallel for collapse(2) schedule(dynamic,1000)
      #pragma omp parallel for schedule(dynamic,1000)
      for (int i = 0; i < vectors->num_vector_local; ++i) {
        for (int fl = 0; fl < vectors->num_field_local; ++fl) {
          vectors_buf->elt<GMFloat>(fl, i) =
            GMVectors_float_get(vectors, fl, i, env);
        }
      }
    } break;
    case MetricType::CCC: {
      // don't use collapse because of overflow for large sizes
      //#pragma omp parallel for collapse(2) schedule(dynamic,1000)
      #pragma omp parallel for schedule(dynamic,1000)
      for (int i = 0; i < vectors->num_vector_local; ++i) {
        for (int fl = 0; fl < vectors->num_packedfield_local; ++fl) {
          vectors_buf->elt<GMBits2x64>(fl, i) =
            GMVectors_bits2x64_get(vectors, fl, i, env);
        }
      }
    } break;
    case MetricType::DUO: {
      // don't use collapse because of overflow for large sizes
      //#pragma omp parallel for collapse(2) schedule(dynamic,1000)
      #pragma omp parallel for schedule(dynamic,1000)
      for (int i = 0; i < vectors->num_vector_local; ++i) {
        for (int fl = 0; fl < vectors->num_packedfield_local; ++fl) {
          vectors_buf->elt<GMBits2x64>(fl, i) =
            GMVectors_bits2x64_get(vectors, fl, i, env);
        }
      }
    } break;
    default:
      COMET_INSIST_INTERFACE(env, false && "Unimplemented metric_type.");
  } // switch
}
//=============================================================================
// checksum of vector entries.

size_t GMVectors_cksum(GMVectors* vectors, CEnv* env) {
  COMET_INSIST(vectors && env);

  if (! env->is_proc_active()) {
    return 0;
  }

  return gm_array_cksum((unsigned char*)(vectors->data), vectors->data_size);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
