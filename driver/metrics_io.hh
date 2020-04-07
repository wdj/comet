//-----------------------------------------------------------------------------
/*!
 * \file   metrics_io.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  I/O utilities for metrics, header.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _comet_metrics_io_hh_
#define _comet_metrics_io_hh_

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

  typedef GMFp32 Float_t;

  template<int N>
  struct Metric {
    uint32_t coords[N];
    Float_t value;
    uint32_t coord0(const CEnv& env) const {
      return env.is_metric_type_bitwise() ? coords[0] / 2 : coords[0];
    }
    uint32_t coord1(const CEnv& env) const {
      return env.is_metric_type_bitwise() ? coords[1] / 2 : coords[1];
    }
    uint32_t coord2(const CEnv& env) const {
      COMET_ASSERT(N >= 3);
      return env.is_metric_type_bitwise() ? coords[2] / 2 : coords[2];
    }
    uint32_t i0(const CEnv& env) const {
      return env.is_metric_type_bitwise() ? coords[0] % 2 : 0;
    }
    uint32_t i1(const CEnv& env) const {
      return env.is_metric_type_bitwise() ? coords[1] % 2 : 0;
    }
    uint32_t i2(const CEnv& env) const {
      COMET_ASSERT(N >= 3);
      return env.is_metric_type_bitwise() ? coords[2] % 2 : 0;
    }
  }; // Metric

  MetricIO(FILE* file, GMMetrics& metrics, CEnv& env);
  ~MetricIO() {}

  void write(size_t coord0, size_t coord1, GMFloat value) const;
  void write(size_t coord0, size_t coord1, size_t coord2,
             GMFloat value) const;
  void write(size_t coord0, size_t coord1, int i0, int i1,
             GMFloat value) const;
  void write(size_t coord0, size_t coord1, size_t coord2,
             int i0, int i1, int i2, GMFloat value) const;

  size_t num_written() const {return num_written_;}

  static size_t num_bytes_written_per_metric(CEnv& env) {
    return env.metric_type() == MetricType::CZEK &&
           env.num_way() == NUM_WAY::_2 ? 4 + 4 + 4 :
           env.metric_type() == MetricType::CZEK &&
           env.num_way() == NUM_WAY::_3 ? 4 + 4 + 4 + 4 :
           env.num_way() == NUM_WAY::_2 ? 4 + 4 + 4 :
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

#endif // _comet_metrics_io_hh_

//-----------------------------------------------------------------------------
