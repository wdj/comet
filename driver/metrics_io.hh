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

//=============================================================================
// Helper class to output metrics values to file.

class MetricIO {
public:

  // NOTE: metrics always written in single precision; this could be changed.

  typedef BasicTypes::FP32 Float_t;

  template<int N>
  struct Metric {
    uint32_t coords[N];
    Float_t value;
    uint32_t iG(const CEnv& env) const {
      return env.is_metric_type_bitwise() ? coords[0] / 2 : coords[0];
    }
    uint32_t jG(const CEnv& env) const {
      return env.is_metric_type_bitwise() ? coords[1] / 2 : coords[1];
    }
    uint32_t kG(const CEnv& env) const {
      COMET_ASSERT(N >= 3);
      return env.is_metric_type_bitwise() ? coords[2] / 2 : coords[2];
    }
    uint32_t iE(const CEnv& env) const {
      return env.is_metric_type_bitwise() ? coords[0] % 2 : 0;
    }
    uint32_t jE(const CEnv& env) const {
      return env.is_metric_type_bitwise() ? coords[1] % 2 : 0;
    }
    uint32_t kE(const CEnv& env) const {
      COMET_ASSERT(N >= 3);
      return env.is_metric_type_bitwise() ? coords[2] % 2 : 0;
    }
  }; // Metric

  MetricIO(FILE* file, GMMetrics& metrics, CEnv& env);
  ~MetricIO() {}

  void write(size_t iG, size_t jG, GMFloat value) const;
  void write(size_t iG, size_t jG, size_t kG, GMFloat value) const;
  void write(size_t iG, size_t jG, int iE, int jE, GMFloat value) const;
  void write(size_t iG, size_t jG, size_t kG,
             int iE, int jE, int kE, GMFloat value) const;

  size_t num_written() const {return num_written_;}

  static size_t num_bytes_written_per_metric(CEnv& env) {
    return env.metric_type() == MetricType::CZEK &&
           env.num_way() == NumWay::_2 ? 4 + 4 + 4 :
           env.metric_type() == MetricType::CZEK &&
           env.num_way() == NumWay::_3 ? 4 + 4 + 4 + 4 :
           env.num_way() == NumWay::_2 ? 4 + 4 + 4 :
                                         4 + 4 + 4 + 4;
  }

  size_t num_bytes_written_per_metric() const {
    return num_bytes_written_per_metric(env_);
  }

  template<int N>
  static void read(Metric<N>& metric, FILE* file, CEnv& env) {

    const size_t num_read = fread(&metric, sizeof(metric), 1, file);
    COMET_INSIST(1 == num_read);
  }

private:

  CEnv& env_;
  FILE* file_;
  const int data_type_;
  int num_way_;
  size_t mutable num_written_;

  template<typename T>
  bool write_(T v, size_t& bytes_written) const {
    const size_t num_written_this = fwrite(&v, sizeof(v), 1, file_);
    bytes_written += sizeof(v);
    return 1 == num_written_this;
  }

  //---Disallowed methods.

  MetricIO(  const MetricIO&);
  void operator=(const MetricIO&);

};

//=============================================================================
// Class to write metrics.

class MetricsIO {

  typedef MetricIO::Float_t Float_t;

public:

  MetricsIO(const char* path_stub, int verbosity, CEnv& env);
  ~MetricsIO();

  void write(GMMetrics& metrics);
  void check_file(GMMetrics& metrics);

  size_t num_written() const {return num_written_;}

  static FILE* open(const char* path_stub, CEnv& env, const char* mode = "wb");

private:

  CEnv& env_;
  const std::string path_stub_;
  FILE* file_;
  int verbosity_;
  size_t num_written_;
  size_t num_written_last_write_;

  size_t bytes_(size_t num) const {
    return num * MetricIO::num_bytes_written_per_metric(env_);
  }

  // Disallowed methods.

  MetricsIO(   const MetricsIO&);
  void operator=(const MetricsIO&);
};

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_METRICS_IO_HH_

//-----------------------------------------------------------------------------
