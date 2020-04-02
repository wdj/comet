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
#define _comet_input_output_hh_

#include "env.hh"
#include "metrics.hh"

#include "driver.hh"

//=============================================================================

namespace comet {

//=============================================================================
// Class to write metrics.

class MetricsIO {
public:

  MetricsIO(const char* path_stub, int verbosity, CEnv& env);
  ~MetricsIO();

  void write(GMMetrics& metrics);

  size_t num_written() {return num_written_;}

  static FILE* open(const char* path_stub, CEnv& env, const char* mode = "w");

private:

  CEnv& env_;
  FILE* file_;
  int verbosity_;
  size_t num_written_;

  // Disallowed methods.

  MetricsIO(   const MetricsIO&);
  void operator=(const MetricsIO&);
};

//=============================================================================
// Helper class to output metrics values to file.

class MetricWriter {
public:

  MetricWriter(FILE* file, GMMetrics* metrics, CEnv* env);
  ~MetricWriter() {}

  size_t get_num_written() {return this->num_written_total_;}
  void write(size_t coord0, size_t coord1, GMFloat value);
  void write(size_t coord0, size_t coord1, size_t coord2, GMFloat value);
  void write(size_t coord0, size_t coord1, int i0, int i1, GMFloat value);
  void write(size_t coord0, size_t coord1, size_t coord2,
             int i0, int i1, int i2, GMFloat value);

private:

  FILE* file_;
  const int data_type_;
  int num_way_;
  CEnv* env_;

  size_t num_written_total_;

  //---Disallowed methods.

  MetricWriter(  const MetricWriter&);
  void operator=(const MetricWriter&);

};

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _comet_metrics_io_hh_

//-----------------------------------------------------------------------------
