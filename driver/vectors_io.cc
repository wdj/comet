//-----------------------------------------------------------------------------
/*!
 * \file   vectors_io.cc
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  I/O utilities for vectors.
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

#include "cstdio"
#include "cstdint"
#include "string.h"

#include "env.hh"
#include "vectors.hh"
#include "vectors_io.hh"

//=============================================================================

namespace comet {

//=============================================================================
/// \brief Read vectors from files: floating point case.

template<typename IO_t = GMFloat>
void static VectorsIO_read_float(GMVectors& vectors, const char* path,
  CEnv& env) {
  COMET_INSIST(path);

  if (! env.is_proc_active())
    return;

  FILE* file = fopen(path, "rb");
  COMET_INSIST(NULL != file && "Unable to open file.");

  const size_t nva = vectors.dm->num_vector_active;
  const size_t nfa = vectors.dm->num_field_active;
  const size_t nvl = vectors.num_vector_local;
  const size_t nfl = vectors.num_field_local;
  const size_t nfal = vectors.dm->num_field_active_local;

  const int fl_min = 0;
  const size_t f_min = fl_min + vectors.dm->field_base;

  // Loop to input vectors

  for (size_t vl = 0; vl < nvl; ++vl) {

    // (global) vector number.
    const size_t v = vl + vectors.dm->vector_base;

    // Fill pad vectors at end with copies of last active vector.
    const size_t v_file = utils::min(v, nva-1);

    // Offset into file of first byte to read
    const size_t loc_min_file = sizeof(IO_t) * (f_min + nfa * v_file);

    // Position to that location in file.
    const int success = fseek(file, loc_min_file, SEEK_SET);
    COMET_INSIST(0 == success && "File seek failure.");

    // Get first location in memory to store into.
    IO_t* const loc_min_mem = GMVectors_float_ptr(&vectors, fl_min, vl, &env);

    // Sanity check on vectors layout in memory.
    // NOTE: the following call is ok since has no side effects
    COMET_INSIST((fl_min+1 >= nfl ||
        GMVectors_float_ptr(&vectors, fl_min+1, vl, &env) == loc_min_mem + 1)
        && "Vector layout is incompatible with operation.");

    const size_t num_to_read = nfal;

    // Perform the read.
    const size_t num_read = fread(loc_min_mem, sizeof(IO_t), num_to_read, file);
    COMET_INSIST(num_read == nfal && "File read failure.");
  } //---for vl

  // Ensure any end pad of each vector is set to zero
  GMVectors_initialize_pad(&vectors, &env);

  fclose(file);
}

//-----------------------------------------------------------------------------
/// \brief Read vectors from files: bits2 case.

void VectorsIO_read_bits2(GMVectors& vectors, const char* path, CEnv& env) {
  COMET_INSIST(path);

  typedef unsigned char IO_t;

  if (! env.is_proc_active())
    return;

  if (0 == vectors.dm->num_field_active_local)
    return;

  FILE* file = fopen(path, "rb");
  COMET_INSIST(NULL != file && "Unable to open file.");

  const size_t nva = vectors.dm->num_vector_active;
  const size_t nfa = vectors.dm->num_field_active;
  const size_t nvl = vectors.num_vector_local;
  const size_t npfl = vectors.dm->num_packedfield_local;
  const size_t nfal = vectors.dm->num_field_active_local;
  const size_t bit_per_f = vectors.dm->num_bit_per_field; // = 2

  const size_t bit_per_pf = vectors.dm->num_bit_per_packedfield;
  const size_t bit_per_byte = 8;
  const size_t f_per_byte = bit_per_byte / bit_per_f;
  const size_t byte_per_pf = bit_per_pf / bit_per_byte;

  const size_t byte_per_v_mem = npfl * byte_per_pf;

  // Stored size of a single vector in the file.
  // NOTE: file format assumes each vector padded up to whole number of bytes.

  const size_t byte_per_v_file = utils::ceil(nfa, f_per_byte);

  const size_t fl_min = 0;
  const size_t f_min = fl_min + vectors.dm->field_base;
  const size_t f_max = f_min + nfal;

  // Buffer to capture read-in data that is possibly not correctly bit-aligned.
  const size_t buf_size = byte_per_v_mem + 1;
  IO_t* const buf = (IO_t*)malloc(buf_size * sizeof(*buf));
  COMET_INSIST(sizeof(IO_t) == 1);

  // Loop to input vectors

  for (size_t vl = 0; vl < nvl; ++vl) {

    // (global) vector number.
    const size_t v = vl + vectors.dm->vector_base;

    // Fill pad vectors at end with copies of last active vector.
    const size_t v_file = utils::min(v, nva-1);

    // Byte in file where this vector starts
    const size_t v_base = byte_per_v_file * v_file;

    // Offset into file to first byte to read
    // NOTE: round downward to capture all bits of required byte.
    const size_t loc_min_file = utils::trunc(f_min, f_per_byte) + v_base;

    // Last field of this vector that must be read (0-based).
    const size_t last_field_needed = f_max - 1;
    COMET_ASSERT(f_max >= 1);

    // Offset into file of last byte of vector to read (0-based).
    const size_t last_byte_needed = utils::trunc(last_field_needed, f_per_byte);

    // Offset into file to [1 plus] last byte to read.
    // NOTE: Rounded upward to capture all bits of required byte.
    const size_t loc_max_file = last_byte_needed + 1 + v_base;
    COMET_INSIST(loc_max_file > loc_min_file);
    COMET_INSIST(loc_max_file <= nva * nfa);

    // Calculate another way, to check.
    COMET_INSIST(loc_max_file == utils::ceil(f_max, f_per_byte) + v_base);

    // Total num bytes to read
    const size_t num_to_read = loc_max_file - loc_min_file;

    // Adjustment factor for first field not falling at beginning of char.
    const int f_offset = f_min % f_per_byte;

    // More checks.
    // NOTE: num bytes read could have at most one additional byte, if straddle.
    COMET_INSIST(num_to_read <= byte_per_v_mem + (f_offset ? 1 : 0));
    COMET_INSIST(num_to_read <= buf_size);

    // Position to that location in file.
    const int success = fseek(file, loc_min_file, SEEK_SET);
    COMET_INSIST(0 == success && "File seek failure.");

    // Perform the read.
    const size_t num_read = fread(buf, sizeof(IO_t), num_to_read, file);
    COMET_INSIST(num_to_read == num_read && "File read failure.");

    // Packed fields object for storing data into vectors struct.
    GMBits2x64 outval = GMBits2x64_null();

    const size_t i_end = utils::min(num_to_read, byte_per_v_mem);

    // Loop over bytes in buffer.

    for (size_t i_in = 0; i_in < i_end; ++i_in) {
      // Get current byte in buf.
      COMET_INSIST(i_in < buf_size);
      const int vlo = buf[i_in];
      // Get next byte in buf.
      const int vhi = (i_in+1 < num_to_read) ? buf[i_in+1] : 0;
      // Create the true byte that is de-offsetted.
      const IO_t inval = ( (vhi << (8-2*f_offset)) |
                           (vlo >> (  2*f_offset)) ) & (int)255;
      // The byte num in the vec to store = the current (lo) byte num of buf.
      const size_t i_out = i_in;
      // Get word of outval to store in.
      const int wordnum = (i_out % 16) / 8;
      // Put into outval.
      outval.data[wordnum] += ((uint64_t)inval) << ((i_out % 8) * 8) ;
      // If outval is full, or now at end of buf, flush outval.
      const bool is_outval_full = i_out % 16 == 15;
      const bool is_end_of_buf = i_in == i_end - 1;
      if (is_outval_full || is_end_of_buf) {
        // Flush buffer
        const size_t pfl = i_out / 16;
        GMVectors_bits2x64_set(&vectors, pfl, vl, outval, &env);
        // Reinitialize outval.
        outval = GMBits2x64_null();
      }
    } // for i

  } // vl

  // Ensure any end of vector pad set to zero
  GMVectors_initialize_pad(&vectors, &env);

  free(buf);
  fclose(file);
}

//-----------------------------------------------------------------------------
/// \brief Read vectors from files.

void VectorsIO::read(GMVectors& vectors, const char* path, CEnv& env) {
  COMET_INSIST(path);

  switch (env.data_type_vectors()) {

    case GM_DATA_TYPE_FLOAT:
      VectorsIO_read_float(vectors, path, env);
    break;

    case GM_DATA_TYPE_BITS2:
      VectorsIO_read_bits2(vectors, path, env);
    break;

    default:
      COMET_INSIST(false && "Invalid data type.");
  } // switch
}

//=============================================================================
/// \brief Write vectors to files: floating point case.

template<typename IO_t = GMFloat>
void static VectorsIO_write_float(GMVectors& vectors, const char* path,
  CEnv& env) {
  COMET_INSIST(path);

  if (! env.is_proc_active())
    return;

  COMET_INSIST_INTERFACE(&env, env.num_proc() == 1 &&
    "Only single proc case currently implemented.");

  FILE* file = fopen(path, "wb");
  COMET_INSIST(NULL != file && "Unable to open file.");

  const size_t nval = vectors.dm->num_vector_active_local;
  const size_t nfal = vectors.dm->num_field_active_local;

  size_t num_written = 0;
  const size_t num_written_attempted = nval * nfal;

  for (size_t vl = 0 ; vl < nval; ++vl) {
    for (size_t fl = 0 ; fl < nfal; ++fl) {

      const IO_t outv = GMVectors_float_get(&vectors, fl, vl, &env);
      const size_t num_to_write = 1;
      const size_t num_written_this = fwrite(&outv, sizeof(IO_t), num_to_write,
                                             file);
      num_written += num_written_this;
    }
  }

  COMET_INSIST(num_written_attempted == num_written &&
    "File write failure.");

  fclose(file);
}

//-----------------------------------------------------------------------------
/// \brief Write vectors to files: bits2 case.

void static VectorsIO_write_bits2(GMVectors& vectors, const char* path,
  CEnv& env) {
  COMET_INSIST(path);

  typedef unsigned char IO_t;

  if (! env.is_proc_active())
    return;

  COMET_INSIST_INTERFACE(&env, env.num_proc() == 1 &&
    "Only single proc case currently implemented.");

  FILE* file = fopen(path, "wb");
  COMET_INSIST(NULL != file && "Unable to open file.");

  const size_t nval = vectors.dm->num_vector_active_local;
  const size_t npfl = vectors.dm->num_packedfield_local;
  const size_t nfal = vectors.dm->num_field_active_local;
  const size_t bit_per_f = vectors.dm->num_bit_per_field; // = 2

  const size_t bit_per_pf = vectors.dm->num_bit_per_packedfield;
  const size_t bit_per_byte = 8;
  const size_t f_per_byte = bit_per_byte / bit_per_f;
  const size_t byte_per_pf = bit_per_pf / bit_per_byte;

  size_t num_written = 0;
  size_t num_written_attempted = 0;

  for (size_t vl = 0 ; vl < nval; ++vl) {
    for (size_t pfl = 0 ; pfl < npfl; ++pfl) {

      GMBits2x64 val = GMVectors_bits2x64_get(&vectors, pfl, vl, &env);

      // Calculate begin and end byte numbers of vector to write

      const size_t offset_min = pfl * byte_per_pf;

      const size_t byte_max = utils::ceil(nfal, f_per_byte);

      const size_t offset_max = utils::min((pfl+1) * byte_per_pf, byte_max);

      // Loop over bytes to output

      for (size_t offset = offset_min, i = 0; offset < offset_max;
           ++offset, ++i) {

        const int index0 = i % 8;
        const int index1 = i / 8;

        const IO_t outval = ((val.data[index1] >> (8*index0)) &
                                (uint64_t)255);

        const size_t num_to_write = 1;

        const size_t num_written_this = fwrite(&outval, sizeof(IO_t),
                                               num_to_write, file);

        num_written += num_written_this;
        num_written_attempted += num_to_write;
      }

    } // pfl
  } // vl

  COMET_INSIST(num_written_attempted == num_written &&
    "File write failure.");

  fclose(file);
}

//-----------------------------------------------------------------------------
/// \brief Write vectors to files.

void VectorsIO::write(GMVectors& vectors, const char* path, CEnv& env) {
  COMET_INSIST(path);

  switch (env.data_type_vectors()) {

    case GM_DATA_TYPE_FLOAT:
      VectorsIO_write_float(vectors, path, env);
    break;

    case GM_DATA_TYPE_BITS2:
      VectorsIO_write_bits2(vectors, path, env);
    break;

    default:
      COMET_INSIST(false && "Invalid data type.");
  } // switch
}

//=============================================================================

void VectorsIO::print(GMVectors& vectors, CEnv& env) {

  if (! env.is_proc_active())
    return;

  const int nval = vectors.dm->num_vector_active_local;
  const int nfal = vectors.dm->num_field_active_local;

  switch (env.data_type_vectors()) {

    case GM_DATA_TYPE_FLOAT: {

      for (int vl = 0; vl < nval; ++vl) {
        for (int fl = 0; fl < nfal; ++fl) {
          const GMFloat value = GMVectors_float_get(&vectors, fl, vl, &env);
            printf("vec_proc %i vec %i field_proc %i field %i value %.16e\n",
                   env.proc_num_vector(), vl,
                   env.proc_num_field(), fl, (double)value);
        } // fl
      } // vl
    } break;

    case GM_DATA_TYPE_BITS2: {

      for (int vl = 0; vl < nval; ++vl) {
        for (int fl = 0; fl < nfal; ++fl) {
          const GMBits2 value = GMVectors_bits2_get(&vectors, fl, vl, &env);
            printf("vec_proc %i vec %i "
                   "field_proc %i field %i value %.1i%.1i\n",
                   env.proc_num_vector(), vl,
                   env.proc_num_field(), fl, value / 2, value % 2);
        } // fl
      }   // vl
    } break;

    default:

      COMET_INSIST(false && "Invalid data_type_vectors.");

  } // switch
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
