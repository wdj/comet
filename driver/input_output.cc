//-----------------------------------------------------------------------------
/*!
 * \file   input_output.cc
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  I/O functions used by driver.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#include "stdio.h"
//#include "stdlib.h"
//#include "stddef.h"
#include "string.h"
//#include "float.h"
//#include "errno.h"

#include "env.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "driver.hh"

//=============================================================================
// Input vectors from files

void set_vectors_from_file(GMVectors* vectors, DriverOptions* do_, GMEnv* env) {
  GMInsist(vectors && do_ && env);
  GMInsist(do_->input_file_path);

  if (! GMEnv_is_proc_active(env)) {
    return;
  }

  switch (GMEnv_data_type_vectors(env)) {
    /*--------------------*/
    case GM_DATA_TYPE_FLOAT: {
    /*--------------------*/

      const int fl = 0;
      const size_t field_base = fl +
        vectors->num_field_local * (size_t)GMEnv_proc_num_field(env);

      FILE* input_file = fopen(do_->input_file_path, "r");
      GMInsist(NULL != input_file && "Unable to open input file.");

      for (int vl = 0; vl < vectors->num_vector_local; ++vl) {

        const size_t proc_num = GMEnv_proc_num_vector_i(env);
        const size_t vector = vl + vectors->num_vector_local * proc_num;
        //---shuffle.
        //const size_t vector = proc_num + GMEnv_num_proc_vector_i(env) * vl;
        /*---Fill pad vectors with copies of the last vector---*/
        const size_t vector_capped = vector <= do_->num_vector_active-1 ?
                                     vector : do_->num_vector_active-1;

        const size_t elt_num = field_base + vectors->num_field * vector_capped;
        const size_t addr_file = elt_num * sizeof(GMFloat);
        int fseek_success = fseek(input_file, addr_file, SEEK_SET);
        fseek_success += 0; /*---Avoid unused var warning---*/
        GMInsist(0 == fseek_success);
        GMFloat* const addr_mem = GMVectors_float_ptr(vectors, fl, vl, env);
        /*---NOTE: the following call is ok since has no side effects---*/
        GMInsist((fl+1 >= vectors->num_field_local ||
            GMVectors_float_ptr(vectors, fl+1, vl, env) == addr_mem + 1)
            && "Vector layout is incompatible with operation.");

        size_t num_read = fread(addr_mem, sizeof(GMFloat),
                                vectors->num_field_local, input_file);
        num_read += 0; /*---Avoid unused var warning---*/
        GMInsist((size_t)vectors->num_field_local == (size_t)num_read);

      } /*---vl---*/

      fclose(input_file);

    } break;
    /*--------------------*/
    case GM_DATA_TYPE_BITS2: {
    /*--------------------*/

      GMInsistInterface(env, GMEnv_num_proc_field(env) == 1 &&
                        "CCC file read for this case not yet implemented.");
      const int pvfl = 0;

      typedef char input_t;

      FILE* input_file = fopen(do_->input_file_path, "r");
      GMInsist(NULL != input_file && "Unable to open input file.");

      for (int vl = 0; vl < vectors->num_vector_local; ++vl) {

        const size_t proc_num = GMEnv_proc_num_vector_i(env);
        const size_t vector = vl + vectors->num_vector_local * proc_num;
        /*---Fill pad vectors with copies of the last vector---*/
        const size_t vector_capped = vector <= do_->num_vector_active-1 ?
                                     vector : do_->num_vector_active-1;

        const int bits_per_field = 2;
        const int bits_per_byte = 8;
	const size_t bytes_per_vector
          = gm_ceil_i8(vectors->num_field * bits_per_field, bits_per_byte);
        const size_t addr_file = bytes_per_vector * vector_capped;

        int fseek_success = fseek(input_file, addr_file, SEEK_SET);
        fseek_success += 0; /*---Avoid unused var warning---*/
        GMInsist(0 == fseek_success);

        input_t* const addr_mem
           = (input_t*)GMVectors_bits2x64_ptr(vectors, pvfl, vl, env);

        size_t num_read = fread(addr_mem, sizeof(input_t),
                                bytes_per_vector, input_file);
        num_read += 0; /*---Avoid unused var warning---*/
        GMInsist(bytes_per_vector == (size_t)num_read);

      } /*---vl---*/

      fclose(input_file);

      //GMInsistInterface(env, false && "Not yet implemented.");
    } break;
    /*--------------------*/
    default:
    /*--------------------*/
      GMInsist(false && "Invalid data type.");
  } /*---switch---*/
}

//=============================================================================
// Class to help output the result metrics values to file

class GMWriter {

  FILE* file_;
  const int data_type_;
  int num_way_;
  GMEnv* env_;

  size_t num_written_total_;

  //typedef unsigned char out_t;
  //const GMFloat out_max = 255;
  //const int buf_size = 4096;
  //int buf_elts = 0;
  //out_t buf[buf_size];

public:

  //--------------------

  GMWriter(FILE* file, GMMetrics* metrics, GMEnv* env) :
    file_(file),
    data_type_(GMEnv_data_type_metrics(env)),
    num_way_(GMEnv_num_way(env)) {

    if (file_ != stdout) {
      GMInsistInterface(env, metrics->num_vector_active ==
                    (GMUInt32)metrics->num_vector_active &&
                    "Too many vectors for output format.");
    }

  }

  //--------------------

  ~GMWriter() {

      //if (file != stdout) {
        //size_t num_written = fwrite(&buf, sizeof(out_t), buf_elts, file);
        //num_written_total += num_written;
        //GMAssert(num_written == buf_elts*sizeof(out_t));
        //buf_elts = 0;
        //printf("Wrote %lu elements of %lu from proc %i.\n",
        //       num_written_total, metrics->num_elts_local,
        //       GMEnv_proc_num(env));
      //}

  }

  //--------------------
  // CZEK 2-way

  void write(size_t coord0, size_t coord1, GMFloat value) {

    bool success = true;

    const GMUInt32 outc0 = coord0;
    size_t num_written = fwrite(&outc0, sizeof(outc0), 1, file_);
    success = success && num_written == 1*sizeof(outc0);

    const GMUInt32 outc1 = coord1;
    num_written = fwrite(&outc1, sizeof(outc1), 1, file_);
    success = success && num_written == 1*sizeof(outc1);

    const GMFp32 outv = value;
    num_written = fwrite(&outv, sizeof(outv), 1, file_);
    success = success && num_written == 1*sizeof(outv);

    num_written_total_ += success ? 1 : 0;

    //out_t out_v = (out_t)(value * out_max);
    //buf[buf_elts++] = out_v;
    //if (buf_elts == buf_size) {
    //  size_t num_written = fwrite(&buf, sizeof(out_t),
    //                              buf_elts, file);
    //  num_written_total += num_written;
    //  GMAssert(num_written == buf_elts*sizeof(out_t));
    //  buf_elts = 0;
    //}
  }

  //--------------------
  // CZEK 3-way

  void write(size_t coord0, size_t coord1, size_t coord2, GMFloat value) {

    bool success = true;

    const GMUInt32 outc0 = coord0;
    size_t num_written = fwrite(&outc0, sizeof(outc0), 1, file_);
    success = success && num_written == 1*sizeof(outc0);

    const GMUInt32 outc1 = coord1;
    num_written = fwrite(&outc1, sizeof(outc1), 1, file_);
    success = success && num_written == 1*sizeof(outc1);

    const GMUInt32 outc2 = coord2;
    num_written = fwrite(&outc2, sizeof(outc2), 1, file_);
    success = success && num_written == 1*sizeof(outc2);

    const GMFp32 outv = value;
    num_written = fwrite(&outv, sizeof(outv), 1, file_);
    success = success && num_written == 1*sizeof(outv);

    num_written_total_ += success ? 1 : 0;
  }

  //--------------------
  // CCC 2-way

  void write(size_t coord0, size_t coord1, int i0, int i1, GMFloat value) {

    bool success = true;

    const GMUInt32 outc0 = i0 + 2 * coord0;
    size_t num_written = fwrite(&outc0, sizeof(outc0), 1, file_);
    success = success && num_written == 1*sizeof(outc0);

    const GMUInt32 outc1 = i1 + 2 * coord1;
    num_written = fwrite(&outc1, sizeof(outc1), 1, file_);
    success = success && num_written == 1*sizeof(outc1);

    const GMFp32 outv = value;
    num_written = fwrite(&outv, sizeof(outv), 1, file_);
    success = success && num_written == 1*sizeof(outv);

    num_written_total_ += success ? 1 : 0;

  }

  //--------------------
  // CCC 3-way

  void write(size_t coord0, size_t coord1, size_t coord2,
             int i0, int i1, int i2, GMFloat value) {

    bool success = true;

    const GMUInt32 outc0 = i0 + 2 * coord0;
    size_t num_written = fwrite(&outc0, sizeof(outc0), 1, file_);
    success = success && num_written == 1*sizeof(outc0);

    const GMUInt32 outc1 = i1 + 2 * coord1;
    num_written = fwrite(&outc1, sizeof(outc1), 1, file_);
    success = success && num_written == 1*sizeof(outc1);

    const GMUInt32 outc2 = i2 + 2 * coord2;
    num_written = fwrite(&outc2, sizeof(outc2), 1, file_);
    success = success && num_written == 1*sizeof(outc2);

    const GMFp32 outv = value;
    num_written = fwrite(&outv, sizeof(outv), 1, file_);
    success = success && num_written == 1*sizeof(outv);

    num_written_total_ += success ? 1 : 0;

  }
};

//=============================================================================
// Output results metrics to file: implementation

void output_metrics_impl(GMMetrics* metrics, DriverOptions* do_,
                         FILE* file, double threshold, GMEnv* env) {
  GMInsist(metrics && do_ && file && env);

  if (! GMEnv_is_proc_active(env)) {
    return;
  }

  /*---Due to redundancy, only results from some processors are needed---*/
  if (GMEnv_proc_num_field(env) != 0) {
    return;
  }

  GMWriter writer(file, metrics, env);

  switch (GMEnv_data_type_metrics(env)) {
    /*--------------------*/
    case GM_DATA_TYPE_FLOAT: {
    /*--------------------*/

      /*----------*/
      if (GMEnv_num_way(env) == GM_NUM_WAY_2) {
      /*----------*/
        size_t index = 0;
        for (index = 0; index < metrics->num_elts_local; ++index) {
          const size_t coord0 =
            GMMetrics_coord_global_from_index(metrics, index, 0, env);
          const size_t coord1 =
            GMMetrics_coord_global_from_index(metrics, index, 1, env);
          if (coord0 >= metrics->num_vector_active ||
              coord1 >= metrics->num_vector_active) {
            continue;
          }
          const GMFloat value
            = GMMetrics_czek_get_from_index(metrics, index, env);
          if (!(threshold < 0. || value > threshold)) {
            continue;
          }
          /*---Output the value---*/
          if (file == stdout) {

            fprintf(file,
              sizeof(GMFloat) == 8 ?
              "element (%li,%li): value: %.17e\n" :
              "element (%li,%li): value: %.8e\n", coord0, coord1, value);

          } else {

            writer.write(coord0, coord1, value);

          }
        }
      }

      /*----------*/
      if (GMEnv_num_way(env) == GM_NUM_WAY_3) {
      /*----------*/
        size_t index = 0;
        for (index = 0; index < metrics->num_elts_local; ++index) {
          const size_t coord0 =
            GMMetrics_coord_global_from_index(metrics, index, 0, env);
          const size_t coord1 =
            GMMetrics_coord_global_from_index(metrics, index, 1, env);
          const size_t coord2 =
            GMMetrics_coord_global_from_index(metrics, index, 2, env);
          if (coord0 >= metrics->num_vector_active ||
              coord1 >= metrics->num_vector_active ||
              coord2 >= metrics->num_vector_active) {
            continue;
          }
          const GMFloat value
            = GMMetrics_czek_get_from_index(metrics, index, env);
          if (!(threshold < 0. || value > threshold)) {
            continue;
          }
          /*---Output the value---*/
          if (file == stdout) {

            fprintf(file,
              sizeof(GMFloat) == 8 ?
              "element (%li,%li,%li): value: %.17e\n" :
              "element (%li,%li,%li): value: %.8e\n",
              coord0, coord1, coord2, value);

          } else {

            writer.write(coord0, coord1, coord2, value);

          }
        }
      }

    } break;
    /*--------------------*/
    case GM_DATA_TYPE_TALLY2X2: {
    /*--------------------*/

      GMInsist(GMEnv_num_way(env) == GM_NUM_WAY_2);

      if (file != stdout) {

        size_t num_buf_ind = 1000000;
        size_t num_buf = 4 * num_buf_ind;

        bool do_out_buf[num_buf];
        int coord0_buf[num_buf];
        int coord1_buf[num_buf];
        int i01_buf[num_buf];
        //size_t loc_buf[num_buf];
        GMFloat value_buf[num_buf];

        //const size_t nva = metrics->num_vector_active;
        //const size_t no_out = 4 * nva * nva;
  
        for (size_t ind_base = 0; ind_base < metrics->num_elts_local;
             ind_base += num_buf_ind) {
          const size_t ind_max = gm_min_i8(metrics->num_elts_local,
                                           ind_base + num_buf_ind);
#pragma omp parallel for collapse(3) schedule(dynamic,10000)
          for (size_t index = ind_base; index < ind_max; ++index) {
            for (int i0 = 0; i0 < 2; ++i0) {
              for (int i1 = 0; i1 < 2; ++i1) {
                const GMFloat value =
                  GMMetrics_ccc_get_from_index_2(metrics, index, i0, i1, env);
                const bool do_thresh = threshold<0. || value>threshold;
                const size_t coord0 = (!do_thresh) ? 0 :
                  GMMetrics_coord_global_from_index(metrics, index, 0, env);
                const size_t coord1 = (!do_thresh) ? 0 :
                  GMMetrics_coord_global_from_index(metrics, index, 1, env);
                const bool do_out = do_thresh &&
                  coord0 < metrics->num_vector_active &&
                  coord1 < metrics->num_vector_active;
                const size_t ind_buf = i1 + 2*(i0 + 2*(index-ind_base));

                do_out_buf[ind_buf] = do_out;
                coord0_buf[ind_buf] = coord0;
                coord1_buf[ind_buf] = coord1;
                //i0_buf[ind_buf] = i0;
                //i1_buf[ind_buf] = i1;
                i01_buf[ind_buf] = i0 + 2*i1;

                //loc_buf[ind_buf] = (!do_out) ? no_out :
                //                   i0 +2*(i1 + 2*(coord0 + nva*coord1));

                value_buf[ind_buf] = value;
              } /*---i1---*/
            } /*---i0---*/
          } /*---index---*/

          const size_t ind_buf_max = 4 * (ind_max - ind_base);
          for (size_t ind_buf = 0; ind_buf < ind_buf_max; ++ind_buf) {
            const bool do_out = do_out_buf[ind_buf];
            //const bool do_out = i01_buf[ind_buf] != -1;
            //const bool do_out = loc_buf[ind_buf] != no_out;
            if (do_out) {
                const int i0 = i01_buf[ind_buf] % 2;
                const int i1 = i01_buf[ind_buf] / 2;
                const size_t coord0 = coord0_buf[ind_buf];
                const size_t coord1 = coord1_buf[ind_buf];
                //const int i0 = loc_buf[ind_buf] % 2;
                //const int i1 = (loc_buf[ind_buf] / 2 ) % nva;
                //const size_t coord0 = (loc_buf[ind_buf] / 4) % nva;
                //const size_t coord1 =  loc_buf[ind_buf] / (4 * nva);
                writer.write(coord0, coord1, i0, i1, value_buf[ind_buf]);
            }
          } /*---ind_buf---*/
        } /*---ind_base---*/

      } else /*---stdout---*/ {

        size_t index = 0;
        for (index = 0; index < metrics->num_elts_local; ++index) {
          const size_t coord0 =
            GMMetrics_coord_global_from_index(metrics, index, 0, env);
          if (coord0 >= metrics->num_vector_active) {
            continue;
          }
          const size_t coord1 =
            GMMetrics_coord_global_from_index(metrics, index, 1, env);
          if (coord1 >= metrics->num_vector_active) {
            continue;
          }
          int num_out_this_line = 0;
          int i0;
          for (i0 = 0; i0 < 2; ++i0) {
            int i1;
            for (i1 = 0; i1 < 2; ++i1) {
              const GMFloat value
                = GMMetrics_ccc_get_from_index_2(metrics, index, i0, i1, env);
              if (!(threshold < 0. || value > threshold)) {
                continue;
              }

              /*---Output the value---*/
              if (file == stdout) {

                if (num_out_this_line == 0) {
                  fprintf(file,
                    "element (%li,%li): values:", coord0, coord1);
                }

                fprintf(file,
                  sizeof(GMFloat) == 8 ?
                  " %i %i %.17e" :
                  " %i %i %.8e", i0, i1, value);

              } else {

                writer.write(coord0, coord1, i0, i1, value);

              }

              num_out_this_line++;

            } /*---i1---*/
          } /*---i0---*/
          if (file == stdout) {
            if (num_out_this_line > 0) {
              fprintf(file, "\n");
            }
          }
        } /*---index---*/

      } /*---if---*/

    } break;

    /*--------------------*/
    case GM_DATA_TYPE_TALLY4X2: {
    /*--------------------*/

      size_t index = 0;
      for (index = 0; index < metrics->num_elts_local; ++index) {
        const size_t coord0 =
          GMMetrics_coord_global_from_index(metrics, index, 0, env);
        const size_t coord1 =
          GMMetrics_coord_global_from_index(metrics, index, 1, env);
        const size_t coord2 =
          GMMetrics_coord_global_from_index(metrics, index, 2, env);
        if (coord0 >= metrics->num_vector_active ||
            coord1 >= metrics->num_vector_active ||
            coord2 >= metrics->num_vector_active) {
          continue;
        }
        int num_out_this_line = 0;
        int i0;
        for (i0 = 0; i0 < 2; ++i0) {
          int i1;
          for (i1 = 0; i1 < 2; ++i1) {
            int i2;
            for (i2 = 0; i2 < 2; ++i2) {
              const GMFloat value
                = GMMetrics_ccc_get_from_index_3(metrics, index, i0, i1, i2,
                                                env);
              if (!(threshold < 0. || value > threshold)) {
                continue;
              }

              /*---Output the value---*/
              if (file == stdout) {

                if (num_out_this_line == 0) {
                  fprintf(file,
                    "element (%li,%li,%li): values:", coord0, coord1, coord2);
                }

                fprintf(file,
                  sizeof(GMFloat) == 8 ?
                  " %i %i %i %.17e" :
                  " %i %i %i %.8e", i0, i1, i2, value);

              } else {

                writer.write(coord0, coord1, coord2, i0, i1, i2, value);

              }

              num_out_this_line++;

            } /*---i2---*/
          } /*---i1---*/
        } /*---i0---*/
        if (file == stdout) {
          if (num_out_this_line > 0) {
            fprintf(file, "\n");
          }
        }
      } /*---index---*/

    } break;
    /*--------------------*/
    default:
      GMInsist(false && "Invalid data type.");
  } /*---switch---*/
}

//=============================================================================
// Output results metrics to file

void output_metrics(GMMetrics* metrics, DriverOptions* do_, GMEnv* env) {
  GMInsist(metrics && do_ && env);

  char* stub = do_->output_file_path_stub;

  /*---Output to file if requested---*/

  if (NULL != stub && GMEnv_is_proc_active(env) &&
      GMEnv_proc_num_field(env) == 0) {

    /*---Form filename---*/

    size_t len = strlen(stub);
    char* path = (char*)malloc((len+50) * sizeof(char));

    GMInsist(env->num_stage < 1000000);
    GMInsist(env->num_phase < 1000000);
    GMInsist(GMEnv_num_proc(env) < 10000000000);

    if (env->num_stage == 1) {
      if (env->num_phase == 1) {
        sprintf(path, "%s_%010i", stub, GMEnv_proc_num(env));
      } else {
        sprintf(path, "%s_%06i_%010i", stub, env->phase_num,
                GMEnv_proc_num(env));
      }
    } else {
      if (env->num_phase == 1) {
        sprintf(path, "%s_%06i_%010i", stub, env->stage_num,
                GMEnv_proc_num(env));
      } else {
        sprintf(path, "%s_%06i_%06i_%010i", stub, env->phase_num,
                env->stage_num, GMEnv_proc_num(env));
      }
    }

    /*---Do output---*/

    double threshold = do_->threshold;

    FILE* file = fopen(path, "w");
    output_metrics_impl(metrics, do_, file, threshold, env);
    fclose(file);
    free(path);
  }

  /*---Output to stdout if requested---*/

  if (do_->verbosity > 1) {
    double threshold = do_->verbosity > 2 ? -1. : do_->threshold;
    output_metrics_impl(metrics, do_, stdout, threshold, env);
  }
}

//=============================================================================
