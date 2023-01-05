//-----------------------------------------------------------------------------
/*!
 * \file   vectors.hh
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

#ifndef _COMET_VECTORS_HH_
#define _COMET_VECTORS_HH_

#include "env.hh"
#include "decomp_mgr.hh"
#include "mirrored_buf.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/*!
 * \class Vectors
 * \brief Store and manipulate sets of vectors used for computing metrics.
 *
 */
//-----------------------------------------------------------------------------

class Vectors {

public:

  // Constructors, destructors.
  Vectors(GMDecompMgr& dm, CEnv& env);
  ~Vectors();

  void allocate();
  void allocate_with_buf();
  void deallocate();

  void initialize();
  void initialize_pad();

  size_t cksum() const;
  void to_buf(MirroredBuf& vectors_buf) const;

  // Accessor functions.

  int num_vector_local() const {return num_vector_local_;}

  int num_packedfield_local() const {return num_packedfield_local_;}

  size_t num_packedfield_vector_local() const {
    return num_packedfield_vector_local_;
  }

  void* __restrict__ data() const {return data_;}

  bool has_buf() const {return has_buf_;}

  GMDecompMgr* dm() const {return &dm_;}

  MirroredBuf* buf() const {
    COMET_INSIST(has_buf_);
    return buf_;
  }

  // Element accessors.

  template<typename Float_t>
  inline Float_t& elt_float(int packedfield_local, int vector_local) {
    COMET_ASSERT(env_.is_double_prec() == (8 == sizeof(Float_t)));
    return elt_<DataTypeId::FLOAT, Float_t>(packedfield_local, vector_local);
  }

  template<typename Float_t>
  inline Float_t elt_float_const(int packedfield_local,
                                 int vector_local) const {
    COMET_ASSERT(env_.is_double_prec() == (8 == sizeof(Float_t)));
    return elt_const_<DataTypeId::FLOAT, Float_t>(packedfield_local,
                                                  vector_local);
  }

  inline GMBits2x64& elt_bits2x64(int packedfield_local, int vector_local) {
    return elt_<DataTypeId::BITS2, GMBits2x64>(packedfield_local, vector_local);
  }

  inline GMBits2x64 elt_bits2x64_const(int packedfield_local,
                                       int vector_local) const {
    return elt_const_<DataTypeId::BITS2, GMBits2x64>(packedfield_local,
                                                     vector_local);
  }

  __host__ __device__
  static size_t index(int packedfield_local,
                      int vector_local,
                      size_t num_packedfield_local_size_t_) {
    // NOTE this function accesses an entire packed value.
    const size_t index = packedfield_local +
      num_packedfield_local_size_t_ * vector_local;
    return index;
  }

  inline GMBits2 bits2_get(int field_local, int vector_local, CEnv& env) const;

  inline void bits2_set(int field_local, int vector_local, GMBits2 value,
                        CEnv& env);

private:

  int num_field_local_;
  int num_vector_local_;
  int num_packedfield_local_;
  size_t num_packedfield_local_size_t_;
  size_t num_packedfield_vector_local_;
  void* __restrict__ data_;
  size_t data_size_;
  bool has_buf_;
  MirroredBuf* buf_;
  CEnv& env_;
  GMDecompMgr& dm_;
  bool is_allocated_;
  int data_type_id_;

  enum {HAS_BUF_TRUE_ = true,
        HAS_BUF_FALSE_ = false};

  void allocate_impl_(bool has_buf);

  // Helpers for accessor functions.

  template<int DTI, class T>
  inline T& elt_(int packedfield_local, int vector_local) {
    COMET_ASSERT(DTI == data_type_id_);
    COMET_ASSERT(packedfield_local >= 0);
    COMET_ASSERT(packedfield_local < num_packedfield_local_);
    COMET_ASSERT(vector_local >= 0);
    COMET_ASSERT(vector_local < num_vector_local_);

    return elt_<DTI,T>(index(packedfield_local, vector_local,
                             num_packedfield_local_));
  }

  template<int DTI, class T>
  inline T elt_const_(int packedfield_local, int vector_local) const {
    COMET_ASSERT(DTI == data_type_id_);
    COMET_ASSERT(packedfield_local >= 0);
    COMET_ASSERT(packedfield_local < num_packedfield_local_);
    COMET_ASSERT(vector_local >= 0 && vector_local < num_vector_local_);

    return elt_const_<DTI,T>(index(packedfield_local, vector_local,
                             num_packedfield_local_));
  }

  template<int DTI, class T>
  T& elt_(size_t index) {
    COMET_ASSERT(DTI == data_type_id_);
    COMET_ASSERT(index+1 >= 0+1);
    COMET_ASSERT(index <
      num_vector_local_ * static_cast<size_t>(num_packedfield_local_));

    return static_cast<T*>(data_)[index];
  }

  template<int DTI, class T>
  T elt_const_(size_t index) const {
    COMET_ASSERT(DTI == data_type_id_);
    COMET_ASSERT(index+1 >= 0+1);
    COMET_ASSERT(index <
      num_vector_local_ * static_cast<size_t>(num_packedfield_local_));

    return static_cast<T*>(data_)[index];
  }
}; // class Vectors

//-----------------------------------------------------------------------------
/*!
 * \brief Accessor: get a single 2-bit (seminibble) value.
 *
 */
GMBits2 Vectors::bits2_get(int field_local,
                             int vector_local,
                             CEnv& env) const {
  COMET_ASSERT(field_local >= 0 && field_local < num_field_local_);
  COMET_ASSERT(vector_local >= 0 && vector_local < num_vector_local_);
  COMET_ASSERT(env.data_type_vectors() == DataTypeId::BITS2);

  // The field address is expressible as a tuple:
  // 1. which GMBits2x64 value,
  // 2. which of the 2 (size2) vector data entries,
  // 3. which of the 32 (size1) 2-bit (size0) fields in the data entry

  const int size0 = 2;
  const int size1 = 32;
  const int size2 = 2;

  int field_index0 = field_local % size1;
  int field_index1 = (field_local / size1) % size2;
  size_t field_index2 = field_local / (size1 * size2);

  const auto ind = index(field_index2, num_packedfield_local_size_t_,
                         vector_local);
  GMBits1_2x64* const __restrict__ address =
      &((static_cast<GMBits2x64*>(data_))[ind].data[field_index1]);

  const auto b_11 = static_cast<GMBits1_2x64>(3);

  return static_cast<GMBits2>(
    ((*address) >> (size0 * field_index0)) & b_11);
}

//-----------------------------------------------------------------------------
/*!
 * \brief Accessor: set a 2-bit (seminibble) value.
 *
 */
void Vectors::bits2_set(int field_local,
                          int vector_local,
                          GMBits2 value,
                          CEnv& env) {
  // This function sets a single 2-bit value.
  COMET_ASSERT(field_local >= 0 && field_local < num_field_local_);
  COMET_ASSERT(vector_local >= 0 && vector_local < num_vector_local_);
  COMET_ASSERT(value+1 >= 1 && value < (1 << GM_BITS2_MAX_VALUE_BITS));
  COMET_ASSERT(env.data_type_vectors() == DataTypeId::BITS2);

  // The field address is expressible as a tuple:
  // 1. which GMBits2x64 value,
  // 2. which of the 2 (size2) vector data entries,
  // 3. which of the 32 (size1) 2-bit (size0) fields in the data entry

  const int size0 = 2;
  const int size1 = 32;
  const int size2 = 2;

  const int field_index0 = field_local % size1;
  const int field_index1 = (field_local / size1) % size2;
  const size_t field_index2 = field_local / (size1 * size2);

  const auto ind = index(field_index2, num_packedfield_local_size_t_,
                         vector_local);
  GMBits1_2x64* const __restrict__ address =
      &((static_cast<GMBits2x64*>(data_))[ind].data[field_index1]);

  const auto b_11 = static_cast<GMBits1_2x64>(3);
  const auto b_value = static_cast<GMBits1_2x64>(value);

  *address &= ~(b_11 << (size0 * field_index0));
  *address |= b_value << (size0 * field_index0);

  COMET_ASSERT(value == bits2_get(field_local, vector_local, env));
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_VECTORS_HH_

//-----------------------------------------------------------------------------
