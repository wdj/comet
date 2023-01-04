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
/*!
 * \brief Constructor.
 *
 */
Vectors::Vectors(GMDecompMgr& dm, CEnv& env)
  : num_field_local_(0)
  , num_vector_local_(0)
  , num_packedfield_local_(0)
  , num_packedfield_local_size_t_(0)
  , num_packedfield_vector_local_(0)
  , data_(NULL)
  , data_size_(0)
  , has_buf_(false)
  , buf_(NULL)
  , env_(env)
  , dm_(dm) 
  , is_allocated_(false)
  , data_type_id_(env.data_type_vectors()) {
  }

//-----------------------------------------------------------------------------
/*!
 * \brief Destructor.
 *
 */
Vectors::~Vectors() {
  deallocate();
}

//-----------------------------------------------------------------------------
/*!
 * \brief Allocate, case of using mirrored buf.
 *
 */
void Vectors::allocate_with_buf() {
  allocate_impl_(HAS_BUF_TRUE_);
}

//-----------------------------------------------------------------------------
/*!
 * \brief Allocate, case of not using mirrored buf.
 *
 */
void Vectors::allocate() {
  allocate_impl_(HAS_BUF_FALSE_);
}

//-----------------------------------------------------------------------------
/*!
 * \brief Allocate, implementation.
 *
 */
void Vectors::allocate_impl_(bool has_buf) {

  if (!env_.is_proc_active())
    return;

  COMET_INSIST(dm_.is_allocated());
  COMET_INSIST(!is_allocated_);

  has_buf_ = has_buf;
  num_field_local_ = dm_.num_field_local;
  num_vector_local_ = dm_.num_vector_local;

  // Allocation size for vector storage

  num_packedfield_local_ = dm_.num_packedfield_local;
  num_packedfield_local_size_t_ = static_cast<size_t>(num_packedfield_local_);

  num_packedfield_vector_local_ =
      num_packedfield_local_ * dm_.num_vector_local;

  data_size_ = num_packedfield_vector_local_ *
              (dm_.num_bit_per_packedfield / BITS_PER_BYTE);

  // Set up vector storage, mirrored buffer.

  buf_ = new MirroredBuf(env_);

  if (has_buf_) {
    buf_->allocate(num_packedfield_local_, num_vector_local_);
    data_ = buf_->h; // alias vector storage to buf
  } else {
    data_ = utils::malloc(data_size_, env_);
  }

  is_allocated_ = true;

  // Set pad entries to zero

  initialize_pad();
}

//-----------------------------------------------------------------------------
/*!
 * \brief Initialize vector entries to zero.
 *
 */
void Vectors::initialize() {

  if (!env_.is_proc_active())
    return;

  COMET_INSIST(is_allocated_);

  const size_t pfl_min = 0;
  const size_t pfl_max = dm_.num_packedfield_local;

  if (DataTypeId::FLOAT == data_type_id_ && !env_.is_double_prec()) {

    typedef float Float_t;
    const auto zero = static_cast<Float_t>(0);;
    for (int vl = 0; vl < num_vector_local_; ++vl) {
      for (size_t pfl = pfl_min; pfl < pfl_max; ++pfl) {
        elt_float<Float_t>(pfl, vl) = zero;
      }
    }

  } else if (DataTypeId::FLOAT == data_type_id_) { //  && env_.is_double_prec()

    typedef double Float_t;
    const auto zero = static_cast<Float_t>(0);;
    for (int vl = 0; vl < num_vector_local_; ++vl) {
      for (size_t pfl = pfl_min; pfl < pfl_max; ++pfl) {
        elt_float<Float_t>(pfl, vl) = zero;
      }
    }

  } else { // DataTypeId::BITS2 == data_type_id_

    const GMBits2x64 zero = GMBits2x64_null();
    for (int vl = 0; vl < num_vector_local_; ++vl) {
      for (size_t pfl = pfl_min; pfl < pfl_max; ++pfl) {
        elt_bits2x64(pfl, vl) = zero;
      } // for pfl
    }

  } // if (DataTypeId::FLOAT == data_type_id_ && !env_.is_double_prec())
}

//-----------------------------------------------------------------------------
/*!
 * \brief Initialize unused (pad) vector entries to zero.
 *
 * The purpose of this function is to ensure final pad words/bits
 * of each vector are set to zero so that word-wise summations of bits
 * aren't corrupted with bad trailing data
 *
 * TODO: possibly initialize pad vectors as well as pad fields.
 * TODO: possibly check to make sure user has set all vector entries (?).
 *
 */
void Vectors::initialize_pad() {

  if (!env_.is_proc_active())
    return;

  COMET_INSIST(is_allocated_);

  const size_t nfal = dm_.num_field_active_local;

  const size_t pfl_min = nfal / dm_.num_field_per_packedfield;
  const size_t pfl_max = dm_.num_packedfield_local;

  if (DataTypeId::FLOAT == data_type_id_ && !env_.is_double_prec()) {

    typedef float Float_t;
    const auto zero = static_cast<Float_t>(0);;
    for (int vl = 0; vl < num_vector_local_; ++vl) {
      for (size_t pfl = pfl_min; pfl < pfl_max; ++pfl) {
        elt_float<Float_t>(pfl, vl) = zero;
      }
    }

  } else if (DataTypeId::FLOAT == data_type_id_) { //  && env_.is_double_prec()

    typedef double Float_t;
    const auto zero = static_cast<Float_t>(0);;
    for (int vl = 0; vl < num_vector_local_; ++vl) {
      for (size_t pfl = pfl_min; pfl < pfl_max; ++pfl) {
        elt_float<Float_t>(pfl, vl) = zero;
      }
    }

  } else { // DataTypeId::BITS2 == data_type_id_

    const auto zero = GMBits2x64_null();
    const uint64_t allbits = 0xffffffffffffffff;
    for (int vl = 0; vl < num_vector_local_; ++vl) {
      for (size_t pfl = pfl_min; pfl < pfl_max; ++pfl) {
        const size_t fl = 64 * pfl;

        if (fl >= nfal) {

          elt_bits2x64(pfl, vl) = zero;

        } else if (fl + 32 >= nfal) {

          GMBits2x64 val = elt_bits2x64_const(pfl, vl);
          const int shift_dist = 64 - 2*(nfal-fl);
          COMET_ASSERT(shift_dist >= 0 && shift_dist < 64);
          val.data[0] &= allbits >> shift_dist;
          val.data[1] = 0;
          elt_bits2x64(pfl, vl) = val;

        } else if (fl + 64 >= nfal) {

          GMBits2x64 val = elt_bits2x64_const(pfl, vl);
          const int shift_dist = 64 - 2*(nfal-fl-32);
          COMET_ASSERT(shift_dist >= 0 && shift_dist < 64);
          val.data[1] &= allbits >> shift_dist;
          elt_bits2x64(pfl, vl) = val;

        } // if
      } // for pfl
    }

  } // if (DataTypeId::FLOAT == data_type_id_ && !env_.is_double_prec())
}

//-----------------------------------------------------------------------------
/*!
 * \brief Deallocate.
 *
 */
void Vectors::deallocate() {

  if (!is_allocated_)
    return;

  if (!env_.is_proc_active())
    return;

  COMET_INSIST(data_);

  if (!has_buf_)
    utils::free(data_, data_size_, env_);

  delete buf_;
  data_ = NULL;

  is_allocated_ = false;
}

//-----------------------------------------------------------------------------
/*!
 * \brief Copy vectors elements to a specified mirrored buffer.
 *
 */
void Vectors::to_buf(MirroredBuf& vectors_buf) const {

  if (!env_.is_proc_active())
    return;

  if (!env_.is_using_linalg())
    return;

  COMET_INSIST(is_allocated_);

  if (env_.metric_type() == MetricType::CZEK && !env_.is_double_prec()) {

      typedef float Float_t;
      // Here and below don't use collapse because of overflow for large sizes
      //#pragma omp parallel for collapse(2) schedule(dynamic,1000)
      #pragma omp parallel for schedule(dynamic,1000)
      for (int vl = 0; vl < num_vector_local_; ++vl) {
        for (int fl = 0; fl < num_field_local_; ++fl) {
          vectors_buf.elt<Float_t>(fl, vl) = elt_float_const<Float_t>(fl, vl);
        }
      }

  } else if (env_.metric_type() == MetricType::CZEK) {
             // && env_.is_double_prec()

      typedef double Float_t;
      #pragma omp parallel for schedule(dynamic,1000)
      for (int vl = 0; vl < num_vector_local_; ++vl) {
        for (int fl = 0; fl < num_field_local_; ++fl) {
          vectors_buf.elt<Float_t>(fl, vl) = elt_float_const<Float_t>(fl, vl);
        }
      }

  } else { // env_.metric_type() == MetricType::CCC ||
           // env_.metric_type() == MetricType::DUO

      #pragma omp parallel for schedule(dynamic,1000)
      for (int vl = 0; vl < num_vector_local_; ++vl) {
        for (int fl = 0; fl < num_packedfield_local_; ++fl) {
          vectors_buf.elt<GMBits2x64>(fl, vl) = elt_bits2x64_const(fl, vl);
        }
      }

  } // if (env_.metric_type() == MetricType::CZEK)
}

//=============================================================================
/// \brief Compute checksum of vector entries.

size_t Vectors::cksum() const {
  if (!env_.is_proc_active())
    return 0;

  return utils::array_cksum(reinterpret_cast<unsigned char*>(data_),
                            data_size_);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
