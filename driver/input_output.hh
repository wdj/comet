//-----------------------------------------------------------------------------
/*!
 * \file   input_output.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  I/O functions used by driver, header.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _comet_input_output_hh_
#define _comet_input_output_hh_

#include "env.hh"
#include "vectors.hh"
#include "metrics.hh"

#include "driver.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------

void set_vectors_from_file(GMVectors* vectors, DriverOptions* do_, GMEnv* env);

void write_vectors_to_file(GMVectors* vectors, const char* vectors_file_path,
                           GMEnv* env);

//=============================================================================
// Class to help output the result metrics values to file

FILE* gm_metrics_file_open(char* metrics_file_path_stub, GMEnv* env);

class MetricWriter {
public:

  MetricWriter(FILE* file, GMMetrics* metrics, GMEnv* env);

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
  GMEnv* env_;

  size_t num_written_total_;

  //---Disallowed methods.

  MetricWriter(  const MetricWriter&);
  void operator=(const MetricWriter&);

};

//=============================================================================

class MetricsFile {
public:

  MetricsFile(DriverOptions* do_, GMEnv* env);

  ~MetricsFile();

  void write(GMMetrics* metrics, GMEnv* env);

  size_t get_num_written() {return num_written_;}

private:

  FILE* file_;
  int verbosity_;
  double threshold_;
  size_t num_written_;

  //---Disallowed methods.

  MetricsFile(   const MetricsFile&);
  void operator=(const MetricsFile&);
};

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _comet_input_output_hh_

//-----------------------------------------------------------------------------
