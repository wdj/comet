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

//=============================================================================
// Vectors constructor.

GMVectors::GMVectors(CEnv& env)
  : num_field(0)
  , num_field_local(0)
  , num_field_active(0)
  , num_vector(0)
  , num_vector_local(0)
  , num_packedfield_local(0)
  , num_packedfield_vector_local(0)
  , data_type_id(0)
  , pad1(0)
  , data(NULL)
  , data_size(0)
  , has_buf_(false)
  , buf_(NULL)
  , env_(env)
  , dm_(NULL) 
  , is_allocated_(false) {
  }

//-----------------------------------------------------------------------------

void GMVectors::allocate(int data_type_id, GMDecompMgr& dm) {
  allocate_impl_(data_type_id, dm, HAS_BUF_FALSE);
}

//-----------------------------------------------------------------------------

void GMVectors::allocate_with_buf(int data_type_id, GMDecompMgr& dm) {
  allocate_impl_(data_type_id, dm, HAS_BUF_TRUE);
}

//-----------------------------------------------------------------------------

void GMVectors::allocate_impl_(int data_type_id,
                               GMDecompMgr& dm,
                               bool has_buf) {
  if (!env_.is_proc_active())
    return;

  this->data_type_id = data_type_id;
  dm_ = &dm;
  has_buf_ = has_buf;
  num_field = dm_->num_field;
  num_field_active = dm_->num_field_active;
  num_field_local = dm_->num_field_local;
  num_vector = dm_->num_vector;
  num_vector_local = dm_->num_vector_local;

  const int bits_per_byte = 8;

  // Allocation size for vector storage

  num_packedfield_local = dm_->num_packedfield_local;

  num_packedfield_vector_local =
      num_packedfield_local * dm_->num_vector_local;

  data_size = num_packedfield_vector_local *
              (dm_->num_bit_per_packedfield / bits_per_byte);

  // Set up vector storage, mirrored buffer

  buf_ = new MirroredBuf(env_);

  if (has_buf_) {
    buf_->allocate(num_packedfield_local,
                   num_vector_local);
    data = buf_->h; // alias vector storage to buf
  } else {
    data = gm_malloc(data_size, &env_);
  }

  // Set pad entries to zero

  initialize_pad();

  is_allocated_ = true;
}

//=============================================================================
// Set vector entries to zero.

void GMVectors::initialize() {

  if (!env_.is_proc_active())
    return;

  const size_t pfl_min = 0;
  const size_t pfl_max = dm_->num_packedfield_local;

  switch (data_type_id) {
    case DataTypeId::FLOAT: {
      GMFloat zero = 0;
      for (int vl = 0; vl < num_vector_local; ++vl) {
        for (size_t pfl = pfl_min; pfl < pfl_max; ++pfl) {
          GMVectors_float_set(this, pfl, vl, zero, &env_);
        }
      }
    } break;
    case DataTypeId::BITS2: {
      const GMBits2x64 zero = GMBits2x64_null();
      for (int vl = 0; vl < num_vector_local; ++vl) {
        for (size_t pfl = pfl_min; pfl < pfl_max; ++pfl) {
          GMVectors_bits2x64_set(this, pfl, vl, zero, &env_);
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

void GMVectors::initialize_pad() {

  /*---Ensure final pad words/bits of each vector are set to zero so that
       word-wise summations of bits aren't corrupted with bad trailing data---*/

  if (!env_.is_proc_active()) {
    return;
  }

  const size_t nfal = dm_->num_field_active_local;

  const size_t pfl_min = nfal / dm_->num_field_per_packedfield;
  const size_t pfl_max = dm_->num_packedfield_local;

  switch (data_type_id) {
    case DataTypeId::FLOAT: {
      GMFloat zero = 0;
      for (int vl = 0; vl < num_vector_local; ++vl) {
        for (size_t pfl = pfl_min; pfl < pfl_max; ++pfl) {
          GMVectors_float_set(this, pfl, vl, zero, &env_);
        }
      }
    } break;
    case DataTypeId::BITS2: {
      const GMBits2x64 zero = GMBits2x64_null();
      const uint64_t allbits = 0xffffffffffffffff;
      for (int vl = 0; vl < num_vector_local; ++vl) {
        for (size_t pfl = pfl_min; pfl < pfl_max; ++pfl) {
          const size_t fl = 64 * pfl;

          if (fl >= nfal) {

            GMVectors_bits2x64_set(this, pfl, vl, zero, &env_);

          } else if (fl + 32 >= nfal) {

            GMBits2x64 val = GMVectors_bits2x64_get(this, pfl, vl, &env_);
            const int shift_dist = 64 - 2*(nfal-fl);
            COMET_ASSERT(shift_dist >= 0 && shift_dist < 64);
            val.data[0] &= allbits >> shift_dist;
            val.data[1] = 0;
            GMVectors_bits2x64_set(this, pfl, vl, val, &env_);

          } else if (fl + 64 >= nfal) {

            GMBits2x64 val = GMVectors_bits2x64_get(this, pfl, vl, &env_);
            const int shift_dist = 64 - 2*(nfal-fl-32);
            COMET_ASSERT(shift_dist >= 0 && shift_dist < 64);
            val.data[1] &= allbits >> shift_dist;
            GMVectors_bits2x64_set(this, pfl, vl, val, &env_);

          } // if
        } // for pfl
      }
    } break;
    default:
      COMET_INSIST(false && "Invalid vectors data_type_id.");
  } // switch
}

//=============================================================================
// Vectors destructor.

void GMVectors::deallocate() {
  if (!is_allocated_)
    return;

  COMET_INSIST(data || ! env_.is_proc_active());

  if (!env_.is_proc_active())
    return;

  if (!has_buf_) {
    gm_free(data, data_size, &env_);
    data = NULL;
  }

  delete buf_;

  is_allocated_ = false;
}

//-----------------------------------------------------------------------------

GMVectors::~GMVectors() {
  deallocate();
}

//=============================================================================
// Copy vectors to mirrored buffer

void GMVectors::to_buf(MirroredBuf& vectors_buf) const {

  if (!env_.is_using_linalg())
    return;

  // Copy vectors into GPU buffers if needed.

  switch (env_.metric_type()) {
    case MetricType::CZEK: {
      // don't use collapse because of overflow for large sizes
      //#pragma omp parallel for collapse(2) schedule(dynamic,1000)
      #pragma omp parallel for schedule(dynamic,1000)
      for (int i = 0; i < num_vector_local; ++i) {
        for (int fl = 0; fl < num_field_local; ++fl) {
          vectors_buf.elt<GMFloat>(fl, i) =
            GMVectors_float_get(this, fl, i, &env_);
        }
      }
    } break;
    case MetricType::CCC: {
      // don't use collapse because of overflow for large sizes
      //#pragma omp parallel for collapse(2) schedule(dynamic,1000)
      #pragma omp parallel for schedule(dynamic,1000)
      for (int i = 0; i < num_vector_local; ++i) {
        for (int fl = 0; fl < num_packedfield_local; ++fl) {
          vectors_buf.elt<GMBits2x64>(fl, i) =
            GMVectors_bits2x64_get(this, fl, i, &env_);
        }
      }
    } break;
    case MetricType::DUO: {
      // don't use collapse because of overflow for large sizes
      //#pragma omp parallel for collapse(2) schedule(dynamic,1000)
      #pragma omp parallel for schedule(dynamic,1000)
      for (int i = 0; i < num_vector_local; ++i) {
        for (int fl = 0; fl < num_packedfield_local; ++fl) {
          vectors_buf.elt<GMBits2x64>(fl, i) =
            GMVectors_bits2x64_get(this, fl, i, &env_);
        }
      }
    } break;
    default:
      COMET_INSIST_INTERFACE(&env_, false && "Unimplemented metric_type.");
  } // switch
}

//=============================================================================
// checksum of vector entries.

size_t GMVectors::cksum() const {

  if (!env_.is_proc_active())
    return 0;

  return gm_array_cksum((unsigned char*)data, data_size);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
