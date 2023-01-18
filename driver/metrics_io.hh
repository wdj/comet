//-----------------------------------------------------------------------------
/*!
 * \file   metrics_io.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  I/O utilities for metrics, header.
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

#ifndef _COMET_METRICS_IO_HH_
#define _COMET_METRICS_IO_HH_

#include "string"

#include "env.hh"
#include "metrics.hh"

#include "driver.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/*!
 * \class MetricForFileTypes
 * \brief Define key types needed for MetricForFile.
 *
 */
//-----------------------------------------------------------------------------

struct MetricForFileTypes {

  // NOTE: currently metrics are always written in single precision.
  // If need to change, this is the place to do it.
  typedef uint32_t IntForFile_t;
  typedef BasicTypes::FP32 FloatForFile_t;

};

//-----------------------------------------------------------------------------
/*!
 * \class MetricForFile
 * \brief Helper class to hold one metric as stored in file.
 *
 */
//-----------------------------------------------------------------------------

template<int NUM_WAY>
class MetricForFile : public MetricForFileTypes {

  // Storage for the coords and value for the metric as stored in file.
  // NOTE: because of the read call, this class cannot have storage for
  // any other members except these two, in this order.
  IntForFile_t coords_[NUM_WAY];
  FloatForFile_t value_;

  // Helper function.
  bool is_bitwise_(const CEnv& env) const {
    return env.is_metric_type_bitwise();
  }

public:

  MetricForFile()
    : coords_()
    , value_(0)
  {}

  // Accessors.

  NV_t iG(const CEnv& env) const {
    const auto result = is_bitwise_(env) ? coords_[0] / 2 : coords_[0];
    return static_cast<NV_t>(result);
  }
  NV_t jG(const CEnv& env) const {
    const auto result = is_bitwise_(env) ? coords_[1] / 2 : coords_[1];
    return static_cast<NV_t>(result);
  }
  NV_t kG(const CEnv& env) const {
    COMET_ASSERT(NUM_WAY >= NumWay::_3);
    const auto result = is_bitwise_(env) ? coords_[2] / 2 : coords_[2];
    return static_cast<NV_t>(result);
  }
  int iE(const CEnv& env) const {
    const auto result = is_bitwise_(env) ? coords_[0] % 2 : 0;
    return static_cast<int>(result);
  }
  int jE(const CEnv& env) const {
    const auto result = is_bitwise_(env) ? coords_[1] % 2 : 0;
    return static_cast<int>(result);
  }
  int kE(const CEnv& env) const {
    COMET_ASSERT(NUM_WAY >= NumWay::_3);
    const auto result = is_bitwise_(env) ? coords_[2] % 2 : 0;
    return static_cast<int>(result);
  }
  FloatForFile_t value() const {return value_;}

  // Read a single metric from file.
  void read(FILE* file, CEnv& env) {
    const size_t num_read = fread(this, sizeof(*this), 1, file);
    COMET_INSIST(1 == num_read);
}

private:

  // Checks.
  void check_() {
    COMET_STATIC_ASSERT(NUM_WAY == NumWay::_2 || NUM_WAY == NumWay::_3);
    COMET_STATIC_ASSERT(sizeof(MetricForFile<NUM_WAY>) ==
      sizeof(IntForFile_t) * NUM_WAY + sizeof(FloatForFile_t));
  }

  // Disallowed methods.

  MetricForFile(const MetricForFile&);
  void operator=(const MetricForFile&);

}; // MetricForFile

//-----------------------------------------------------------------------------
/*!
 * \class SingleMetricIO
 * \brief Helper class to read/write single metric at a time.
 *
 */
//-----------------------------------------------------------------------------

class SingleMetricIO : public MetricForFileTypes {

public:

  // Constructors, destructors.

  SingleMetricIO(FILE* file, GMMetrics& metrics, CEnv& env);
  ~SingleMetricIO() {}

  // Write a single metric with indexing to file.

  template<typename Float_t>
  void write(NV_t iG, NV_t jG, Float_t value) const;

  template<typename Float_t>
  void write(NV_t iG, NV_t jG, NV_t kG, Float_t value) const;

  template<typename Float_t>
  void write(NV_t iG, NV_t jG, int iE, int jE, Float_t value) const;

  template<typename Float_t>
  void write(NV_t iG, NV_t jG, NV_t kG,
             int iE, int jE, int kE, Float_t value) const;

  // Sizes.

  size_t num_written() const {return num_written_;}

  static size_t num_bytes_written_per_metric(CEnv& env) {
    return env.num_way() * sizeof(IntForFile_t) + sizeof(FloatForFile_t);
  }

  size_t num_bytes_written_per_metric() const {
    return num_bytes_written_per_metric(env_);
  }

private:

  CEnv& env_;
  FILE* file_;
  const int data_type_;
  int num_way_;
  size_t mutable num_written_;

  // Helper function to write a scalar value to the file.

  template<typename T>
  bool write_(T v, size_t& bytes_written) const {
    const size_t num_written_this = fwrite(&v, sizeof(v), 1, file_);
    bytes_written += sizeof(v);
    const bool is_written_correctly = 1 == num_written_this;
    return is_written_correctly;
  }

  friend class MetricsIO;

  // Disallowed methods.

  SingleMetricIO(const SingleMetricIO&);
  void operator=(const SingleMetricIO&);

};

//-----------------------------------------------------------------------------
/*!
 * \class MetricsIO
 * \brief Class to write metrics.
 *
 */
//-----------------------------------------------------------------------------

class MetricsIO : public MetricForFileTypes {

public:

  // Constructors, destructors.

  MetricsIO(const char* path_stub, int verbosity, CEnv& env);
  ~MetricsIO();

  void deallocate();

  void write(GMMetrics& metrics);

  void check_file(GMMetrics& metrics);

  template<typename Float_t>
  void check_file_impl_(GMMetrics& metrics);

  size_t num_written() const {return num_written_;}

  static bool can_write_file(const char* path_stub, CEnv& env);

private:

  // Set this false to close/reopen files between stages/phases.
  bool is_leaving_files_open_() {return true;}

  void open_(const char* mode = "wb");
  void close_();

  CEnv& env_;
  const std::string path_stub_str_;
  FILE* file_;
  int verbosity_;
  NML_t num_written_;
  NML_t num_written_last_write_;
  bool is_allocated_;

  size_t bytes_(size_t num) const {
    return num * SingleMetricIO::num_bytes_written_per_metric(env_);
  }

  bool is_path_stub_() const {return strlen(path_stub_str_.c_str()) > 0;}

  // Disallowed methods.

  MetricsIO(const MetricsIO&);
  void operator=(const MetricsIO&);
};

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_METRICS_IO_HH_

//-----------------------------------------------------------------------------
