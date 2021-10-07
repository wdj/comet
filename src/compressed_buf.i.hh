//-----------------------------------------------------------------------------
/*!
 * \file   compressed_buf.i.hh
 * \author Wayne Joubert
 * \date   Tue May  5 20:11:12 EDT 2020
 * \brief  Mirrored buffer allowing for compression.
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

#ifndef _COMET_COMPRESSED_BUF_I_HH_
#define _COMET_COMPRESSED_BUF_I_HH_

#include <type_traits>

#include "env.hh"
#include "mirrored_buf.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------

template<typename T> T CompressedBufAccessorUnspecialized_<T>::
elt_const(size_t ind0, size_t ind1, const CompressedBuf* cbuf) {
  return cbuf->buf_->elt_const<T>(ind0, ind1);
}

//-----------------------------------------------------------------------------
/// \brief Cartesian elt access, general (uncompressed) case.

template<typename T> T CompressedBufAccessor_<T>::
elt_const(size_t ind0, size_t ind1, const CompressedBuf* cbuf) {
  COMET_ASSERT(!cbuf->do_compress());
  return CompressedBufAccessorUnspecialized_<T>::elt_const(ind0, ind1, cbuf);
  //return cbuf->buf_->elt_const<T>(ind0, ind1);
}

//-----------------------------------------------------------------------------
/// \brief Cartesian elt access, compressed case, single metric entry case.

// NOTE: this needs to be here because Tally2x2<MetricFormat::SINGLE>::TypeIn
// is float, which also is incidentally GMFloat for Czekanowski for
// the single precision build, so that case comes here.  Thus needed to
// properly disambiguate types.

inline Tally2x2<MetricFormat::SINGLE>::TypeIn
CompressedBufAccessor_<Tally2x2<MetricFormat::SINGLE>::TypeIn>::
elt_const(size_t ind0, size_t ind1, const CompressedBuf* cbuf) {
  COMET_ASSERT(cbuf->env_.metric_type() == MetricType::CZEK);
  typedef Tally2x2<MetricFormat::SINGLE>::TypeIn T;
  return CompressedBufAccessorUnspecialized_<T>::elt_const(ind0, ind1, cbuf);
  //return cbuf->buf_->elt_const<T>(ind0, ind1);
  //COMET_INSIST(false && "Case not supported.");
  //typedef Tally2x2<MetricFormat::SINGLE>::TypeIn T;
  //return (T)0;
}

//-----------------------------------------------------------------------------
/// \brief Cartesian elt accessor, compressed case, entire metric (4 entries).

inline Tally2x2<MetricFormat::SINGLE>
CompressedBufAccessor_<Tally2x2<MetricFormat::SINGLE>>::
elt_const(size_t ind0, size_t ind1, const CompressedBuf* cbuf) {
  COMET_ASSERT(CompressedBuf::State::IDLE == cbuf->state_);
  COMET_ASSERT(ind0 < cbuf->buf_->dim0 && ind1 < cbuf->buf_->dim1);

  typedef Tally2x2<MetricFormat::SINGLE> T;
  typedef CompressedBuf CBuf;
  typedef CBuf::MFTTypeIn MFTTypeIn; // = float
  const MirroredBuf* buf_ = cbuf->buf_;
  CBuf::Reader& reader_ = cbuf->reader_;

  // index into Tally2x2<MetricFormat::SINGLE> elt of (uncompresed) buf.
  const size_t ind = ind0 + buf_->dim0 * ind1;

  // NOTE: for compressed case, assuming here that access are in fact sequential.
  COMET_ASSERT(ind > reader_.ind_recent || ! reader_.is_reading_started ||
               !cbuf->do_compress_);

  // Record details of this call.
  reader_.is_reading_started = true;
  reader_.ind_recent = ind;
  reader_.ind0_recent = ind0;
  reader_.ind1_recent = ind1;

  // Trap simple (uncompressed) case.
  //if (!cbuf->do_compress_) {
  //  auto i = buf_->elt_const<T>(ind0, ind1);
  //  printf("%f\n", (double)*(float*)&i); //FIX
  //}
//printf("//////////////   %i %i %i\n", (int)ind0, (int)ind1, (int)buf_->dim0);
  if (!cbuf->do_compress_)
    return buf_->elt_const<T>(ind0, ind1);

  T result = T::null();

# ifdef COMET_ASSERTIONS_ON
    // number of Tally2x2<MetricFormat::SINGLE> elts of (uncompresed) buf.
    const size_t dim = buf_->dim0 * buf_->dim1;
    // number of TypeIn (float) elements of (uncompressed) buf.
    const size_t dim_typein = 4 * dim;
# endif

  // Index to run in the rle.
  size_t& ind_run = reader_.ind_run_recent;

  // Loop to look for and pick up 4 table entries, if they are in rle.

  // NOTE: the order of next 2 nested loops must match mem layout of Tally2x2,
  // since that is the order of elements submitted to the rle.

  for (int iE=0; iE<2; ++iE) {
    for (int jE=0; jE<2; ++jE) {

      // index into TypeIn (float) elts of (uncompresed) buf.
      const size_t ind_typein = jE + 2 * ( iE + 2 * ind );

      COMET_ASSERT(ind_typein >= reader_.lengths_running_sum);
      COMET_ASSERT(ind_typein < dim_typein);

      // Loop over rle, starting at current location, to seek element.

      for( ; ind_run <  cbuf->num_runs_; ) {

          // Get TypeIn (float) elt indices (uncompressed) covered by this run.
          const size_t length_run = cbuf->lengths_alias_buf_.
            elt_const<CBuf::Lengths_t>(ind_run, 0);
          const size_t ind_typein_min = reader_.lengths_running_sum;
          const size_t ind_typein_max = reader_.lengths_running_sum +
                                        length_run;

          if (ind_typein >= ind_typein_min &&
              ind_typein < ind_typein_max) {
            const MFTTypeIn key_run =
              cbuf->keys_alias_buf_.elt_const<MFTTypeIn>(ind_run, 0);
            T::set(result, iE, jE, key_run);
            reader_.ind_typein_recent = ind_typein;
            reader_.iE_recent = iE;
            reader_.jE_recent = jE;
            break;
          }

          ++ind_run;
          reader_.lengths_running_sum += length_run;

      } // for

    } // for jE
  } // for iE

// Turn this on to explicitly check values against uncompressed buf.
//# define _COMET_COMPRESSED_BUF_CHECK_RESULT

# ifdef _COMET_COMPRESSED_BUF_CHECK_RESULT
    for (int iE=0; iE<2; ++iE)
      for (int jE=0; jE<2; ++jE)
        COMET_INSIST(T::get(buf_->elt_const<T>(ind0, ind1), iE, jE) ==
                     T::get(result, iE, jE));
# endif
  //printf("%f\n", (double)*(float*)&result);
  return result;
}

//-----------------------------------------------------------------------------
/// \brief Sequential elt access, general case.

template<typename T> T CompressedBufAccessor_<T>::
elt_const(size_t ind_entry, const CompressedBuf* cbuf) {
  COMET_INSIST(false && "Case not supported.");
  return 0;
}

//-----------------------------------------------------------------------------
/// \brief Sequential elt access, compressed case, entire metric (4 entries).

inline Tally2x2<MetricFormat::SINGLE>
CompressedBufAccessor_<Tally2x2<MetricFormat::SINGLE>>::
elt_const(size_t ind_entry, const CompressedBuf* cbuf) {
  COMET_INSIST(false && "Case not supported.");
  return T::null();;
}

//-----------------------------------------------------------------------------
/// \brief Sequential elt access, compressed case, single metric entry case.

inline Tally2x2<MetricFormat::SINGLE>::TypeIn
CompressedBufAccessor_<Tally2x2<MetricFormat::SINGLE>::TypeIn>::
elt_const(size_t ind_entry, const CompressedBuf* cbuf) {
  COMET_ASSERT(CompressedBuf::State::IDLE == cbuf->state_);
  COMET_ASSERT(ind_entry < cbuf->num_entries_);
  COMET_ASSERT((cbuf->reader_.ind_entry_recent+1 == ind_entry ||
                !cbuf->reader_.is_reading_started ||
                !cbuf->do_compress_) && "Assume sequential access.");

  typedef CompressedBuf CBuf;
  typedef CBuf::TTable_t TTable_t;
  typedef Tally2x2<MetricFormat::SINGLE>::TypeIn TypeIn; // = float
  //typedef CBuf::MFTTypeIn MFTTypeIn; // = float
  const MirroredBuf* buf_ = cbuf->buf_;
  CBuf::Reader& reader_ = cbuf->reader_;

  // Record details of this call.
  reader_.is_reading_started = true;
  reader_.ind_entry_recent = ind_entry;

  size_t& ind_typein = reader_.ind_typein_recent;
  int& iE = reader_.iE_recent;
  int& jE = reader_.jE_recent;
  size_t& ind = reader_.ind_recent;
  size_t& ind0 = reader_.ind0_recent;
  size_t& ind1 = reader_.ind1_recent;

  // Trap simple (uncompressed) case.
  if (!cbuf->do_compress_) {
    ind_typein = ind_entry;
    iE = (ind_typein / 2) % 2;
    jE = ind_typein % 2;
    ind = ind_typein / 4;
    ind0 = ind % buf_->dim0;
    ind1 = ind / buf_->dim0;
    return TTable_t::get(buf_->elt_const<TTable_t>(ind0, ind1), iE, jE);
  }

  size_t& lengths_running_sum = reader_.lengths_running_sum;
  size_t& ind_run = reader_.ind_run_recent;
  size_t& ind_in_run = reader_.ind_in_run_recent;

  // Skip over run of zeros in rle, if present at current location.

  TypeIn key_run = cbuf->keys_alias_buf_.elt_const<TypeIn>(ind_run, 0);

  if (!key_run) {
    const size_t length_run = cbuf->lengths_alias_buf_.
      elt_const<CBuf::Lengths_t>(ind_run, 0);
    lengths_running_sum += length_run;
    ++ind_run;
//printf("%f %i %i\n", (double) key_run, (int)length_run, (int)ind_run);
    COMET_ASSERT(ind_run < cbuf->num_runs_);
    ind_in_run = 0;
    // Now pointing to start of a run of a nonzero value.
    // Select this value to return.
    key_run = cbuf->keys_alias_buf_.elt_const<TypeIn>(ind_run, 0);
    COMET_ASSERT(key_run);
  }

  // Increment to indicate to where to start looking on the next call.

  const size_t length_run = cbuf->lengths_alias_buf_.
    elt_const<CBuf::Lengths_t>(ind_run, 0);
  COMET_ASSERT(ind_in_run < length_run);

  ++ind_in_run;

  if (ind_in_run == length_run) {
    ++ind_run;
    ind_in_run = 0;
  }

  // Get current (uncompressed) location.

  ind_typein = lengths_running_sum;
  ++lengths_running_sum;

  // Finish.

  iE = (ind_typein / 2) % 2;
  jE = ind_typein % 2;
  ind = ind_typein / 4;
  ind0 = ind % buf_->dim0;
  ind1 = ind / buf_->dim0;

//if (! TTable_t::get(buf_->elt_const<TTable_t>(ind0, ind1), iE, jE) == key_run)
//printf("%i %i %i %i %i %i %i  %f  %f\n", (int)ind_entry, (int)ind_typein, iE, jE, (int)ind, (int)ind0, (int)ind1, (double)key_run, (double)TTable_t::get(buf_->elt_const<TTable_t>(ind0, ind1), iE, jE));

  //COMET_ASSERT(TTable_t::get(buf_->elt_const<TTable_t>(ind0, ind1), iE, jE) == key_run);

  return key_run;
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_COMPRESSED_BUF_I_HH_

//-----------------------------------------------------------------------------
