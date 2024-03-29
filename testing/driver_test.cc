//-----------------------------------------------------------------------------
/*!
 * \file   driver_test.cc
 * \author Wayne Joubert
 * \date   Fri Nov  6 18:18:21 EST 2015
 * \brief  Tester for driver.
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

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "string"

#ifdef COMET_USE_GTEST
# include "gtest/gtest.h"
# define BEGIN_TESTS
# define END_TESTS
#else
# define GTEST_API_
# define EXPECT_EQ(a, b) COMET_INSIST((a) == (b));
# define BEGIN_TESTS int RUN_ALL_TESTS() {
# define END_TESTS return 0;}
# define TEST(a,b)
#endif

#include "env.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "checksum.hh"
#include "compute_metrics.hh"

#include "driver.hh"
#include "vectors_io.hh"
#include "metrics_io.hh"
#include "test_problems.hh"

enum {PROCS_MAX = TEST_PROCS_MAX};

//=============================================================================

bool can_run(const char* options) {

  COMET_INSIST(options);

  comet::CEnv env_nocomms(options, PROCS_MAX, 0);

  if (env_nocomms.num_proc() > 1 && ! comet::BuildHas::MPI)
    return false;

  comet::CEnv env(options);

  return env.can_run();

}

//=============================================================================

bool is_using_tc(const char* options) {

  COMET_INSIST(options);

  comet::CEnv env(options);

  return env.is_using_tc();
}

//=============================================================================

bool compare_2runs(const char* options1, const char* options2) {
  COMET_INSIST(options1 && options2);

  if (!(can_run(options1) && can_run(options2)))
    return true;

  // Do 2 runs.

  const bool is_proc_num_0 = comet::System::is_proc_num_0();

  if (is_proc_num_0) {
    printf("%s\n", options1);
  }

  comet::Checksum checksum1;
  comet::Driver::perform_run(checksum1, options1);

  if (is_proc_num_0) {
    printf("%s\n", options2);
  }

  comet::Checksum checksum2;
  comet::Driver::perform_run(checksum2, options2);

  // Need test result only on proc 0.

  const bool is_passed = is_proc_num_0 ? checksum1.is_equal(checksum2) : true;;

  return is_passed;
}

//=============================================================================

bool compare_3runs(const char* options1,
                   const char* options2,
                   const char* options3) {
  COMET_INSIST(options1 && options2 && options3);

  // Test whatever is possible/supported.

  if (!can_run(options1))
    return compare_2runs(options2, options3);

  if (!can_run(options2))
    return compare_2runs(options1, options3);

  if (!can_run(options3))
    return compare_2runs(options1, options2);

  //if (!(can_run(options1) && can_run(options2) && can_run(options3)))
  //  return true;

  // Do 3 runs.

  const bool is_proc_num_0 = comet::System::is_proc_num_0();

  if (is_proc_num_0) {
    printf("%s\n", options1);
  }
  comet::Checksum checksum1;
  comet::Driver::perform_run(checksum1, options1);

  if (is_proc_num_0) {
    printf("%s\n", options2);
  }
  comet::Checksum checksum2;
  comet::Driver::perform_run(checksum2, options2);

  if (is_proc_num_0) {
    printf("%s\n", options3);
  }
  comet::Checksum checksum3;
  comet::Driver::perform_run(checksum3, options3);

  // Need test result only on proc 0.

  const bool is_passed = is_proc_num_0 ?  checksum1.is_equal(checksum2) &&
                                          checksum1.is_equal(checksum3) : true;
  return is_passed;
}

//=============================================================================

void test_2runs(const char* options1,
                const char* options2) {
  EXPECT_EQ(true, compare_2runs(options1, options2));
}

//=============================================================================

void create_vectors_file(const char* file_path, int num_field, int num_vector,
                         int metric_type, int num_way, int problem_type,
                         int verbosity) {

  using namespace comet;

  std::string options = " --num_way " + std::to_string(num_way);
  options += " --metric_type " + std::string(MetricType::str(metric_type));
  options += " --num_proc_vector " + std::to_string(1);
  options += " --compute_method REF";
  options += " --verbosity 1";
  CEnv env_value(MPI_COMM_WORLD, options.c_str());
  CEnv* env = &env_value;

  GMDecompMgr dm(env_value);
  dm.allocate(false, false, num_field, num_vector, env->data_type_vectors());

  Vectors vectors(dm, env_value);
  vectors.allocate();
  TestProblem::set_vectors_synthetic(&vectors, problem_type, verbosity, env);

  VectorsIO::write(vectors, file_path, *env);

  vectors.deallocate();
  dm.deallocate();
}

//=============================================================================

void DriverTest_czek2_() {
#if 1

  //----------
  //---2-way, all2all no
  //----------

#if 0
  EXPECT_EQ(
      true,
      compare_2runs("--metric_type czekanowski --num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method GPU",
                    "--metric_type czekanowski --num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method CPU"));
#endif

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 1 --num_vector_local 2 "
                    "--compute_method CPU",
                    "--num_proc_vector 1 --num_field 1 --num_vector_local 2 "
                    "--compute_method GPU"));

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method CPU",
                    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method GPU"));

  //----------
  //---2-way, all2all yes, small
  //----------

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 1 --num_vector_local 2 "
                    "--compute_method CPU ",
                    "--num_proc_vector 1 --num_field 1 --num_vector_local 2 "
                    "--compute_method CPU --all2all yes"));

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 1 --num_vector_local 2 "
                    "--compute_method CPU ",
                    "--num_proc_vector 1 --num_field 1 --num_vector_local 2 "
                    "--compute_method GPU --all2all yes"));

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 1 --num_vector_local 4 "
                    "--compute_method CPU ",
                    "--num_proc_vector 2 --num_field 1 --num_vector_local 2 "
                    "--compute_method CPU --all2all yes"));

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 1 --num_vector_local 20 "
                    "--compute_method CPU ",
                    "--num_proc_vector 10 --num_field 1 --num_vector_local 2 "
                    "--compute_method CPU --all2all yes"));

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 1 --num_vector_local 4 "
                    "--compute_method CPU ",
                    "--num_proc_vector 2 --num_field 1 --num_vector_local 2 "
                    "--compute_method GPU --all2all yes"));

  EXPECT_EQ(
      true,
      compare_3runs("--num_proc_vector 1 --num_field 2 --num_vector 5 "
                    "--compute_method REF --all2all yes",
                    "--num_proc_vector 2 --num_field 2 --num_vector 5 "
                    "--compute_method CPU --all2all yes",
                    "--num_proc_vector 2 --num_field 2 --num_vector 5 "
                    "--compute_method GPU --all2all yes"));
#endif

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_proc_field 1 "
                    "--num_field 7 --num_vector 5 "
                    "--compute_method REF --all2all yes",
                    "--num_proc_vector 1 --num_proc_field 3 "
                    "--num_field 7 --num_vector 5 "
                    "--compute_method GPU --all2all yes"));

#if 1
  //----------
  //---2-way, all2all yes, large
  //----------

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method CPU",
                    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method CPU --all2all yes"));

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method CPU",
                    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method GPU"
                    " --all2all yes"));

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method CPU",
                    "--num_proc_vector 2 --num_field 100 --num_vector_local 24 "
                    "--compute_method CPU --all2all yes"));

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method CPU",
                    "--num_proc_vector 2 --num_field 100 --num_vector_local 24 "
                    "--compute_method GPU --all2all yes"));

  //----------
  //---num_proc_field
  //----------

  EXPECT_EQ(true, compare_2runs("--num_proc_vector 1 --num_proc_field "
                                "1 --num_field 2 --num_vector_local 2 "
                                "--compute_method CPU",
                                "--num_proc_vector 1 --num_proc_field "
                                "2 --num_field 2 --num_vector_local 2 "
                                "--compute_method GPU"));

  EXPECT_EQ(true, compare_2runs("--num_proc_vector 1 --num_proc_field "
                                "1 --num_field 2 --num_vector_local 4 "
                                "--compute_method CPU",
                                "--num_proc_vector 2 --num_proc_field "
                                "2 --num_field 2 --num_vector_local 2 "
                                "--compute_method GPU --all2all yes"));
#endif

  //----------
  //---num_proc_repl, 2-way
  //----------

  char options1[1024];
  char options2[1024];

  char options_template_1[] =
      "--metric_type czekanowski "
      "--num_field 4 --num_vector_local %i --compute_method %s --all2all yes "
      "--num_proc_vector %i --num_proc_repl %i "
      "--num_proc_field %i --num_way %i --num_stage %i";

  for (int gpu=0; gpu<=1; ++gpu) {
  //for (int gpu=0; gpu<=0; ++gpu) {
    for (int num_vector_local=4; num_vector_local<=5; ++num_vector_local) {
    //for (int num_vector_local=4; num_vector_local<=4; ++num_vector_local) {
      for (int num_proc_vector=1; num_proc_vector<=6; ++num_proc_vector) {
      //for (int num_proc_vector=1; num_proc_vector<=1; ++num_proc_vector) {
        for (int num_proc_repl=1; num_proc_repl<=6; ++num_proc_repl) {
        //for (int num_proc_repl=1; num_proc_repl<=1; ++num_proc_repl) {
          const int num_proc_field = gpu ? 2 : 1;
          //const int num_proc_field = 1;
          if (num_proc_vector * num_proc_field * num_proc_repl > PROCS_MAX) {
            continue;
          }
          const int num_way = 2;
          sprintf(options1, options_template_1,
                  num_vector_local*num_proc_vector,
                  "GPU", 1, 1,
                  1, num_way, 1);
          sprintf(options2, options_template_1, num_vector_local,
                  gpu ? "GPU" : "CPU", num_proc_vector, num_proc_repl,
                  num_proc_field, num_way, 1);
          test_2runs(options1, options2);
        }
      }
    }
  }

  //----------
  //---num_phase, 2-way
  //----------

  for (int num_proc_vector=1; num_proc_vector<=8; ++num_proc_vector) {
    for (int num_proc_repl=1; num_proc_repl<=8; ++num_proc_repl) {
      for (int num_phase=2; num_phase<=8; ++num_phase) {
        if (!(num_phase <= 1 + num_proc_vector/2)) {
          continue;
        }
        char options_template[] =
          "--metric_type czekanowski "
          "--num_field 7 --num_vector 12 --compute_method GPU --all2all yes "
          "--num_proc_vector %i --num_proc_repl %i --num_phase %i "
          "--num_way 2";
        sprintf(options1, options_template, num_proc_vector, num_proc_repl,
                num_phase);
        sprintf(options2, options_template, 1, 1, 1);
        test_2runs(options1, options2);
      }
    }
  }

  //----------
  //---file input
  //----------

  for (int num_vector = 13; num_vector <= 13; ++num_vector) {
    for (int num_field = 1; num_field <= 10; ++num_field) {
      create_vectors_file("czek_2way_in.bin", num_field, num_vector,
                          comet::MetricType::CZEK, 2, comet::ProblemType::DEFAULT, 1);
      for (int num_proc_vector=1; num_proc_vector<=4; ++num_proc_vector) {
        for (int num_proc_field=1; num_proc_field<=5; ++num_proc_field) {

          char options1[1024];
          char options2[1024];

          char options_template[] =
                 "--num_vector %i --num_field %i "
                 "--num_proc_vector %i --num_proc_field %i "
                 "--compute_method GPU --all2all yes %s --verbosity 1";

          sprintf(options1, options_template,
                  num_vector, num_field, num_proc_vector, num_proc_field, "");

          sprintf(options2, options_template,
                  num_vector, num_field, num_proc_vector, num_proc_field,
                  "--input_file czek_2way_in.bin");

          test_2runs(options1, options2);
        }
      }
    }
  }

  //----------
  //---Misc options
  //----------

  // TODO: set up better test
  //EXPECT_EQ(
  //    true,
  //    compare_2runs("--num_proc_vector 1 --num_field 3 --num_vector 3 "
  //                  "--compute_method GPU",
  //                  "--num_proc_vector 1 --num_field 3 --num_vector 3 "
  //                  "--compute_method GPU --checksum no"));
} // DriverTest_czek2_

//=============================================================================

void DriverTest_czek3_() {
#if 1
  //----------
  //---3-way, all2all no
  //----------

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 1 --num_vector_local 3 "
                    "--compute_method CPU "
                    "--num_way 3",
                    "--num_proc_vector 1 --num_field 1 --num_vector_local 3 "
                    "--compute_method GPU --num_way 3"));
  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method CPU --num_way 3",
                    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method GPU --num_way 3"));

  //----------
  //---3-way, all2all yes, small
  //----------

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 1 --num_vector_local 3 "
                    "--compute_method CPU --num_way 3",
                    "--num_proc_vector 1 --num_field 1 --num_vector_local 3 "
                    "--compute_method CPU --num_way 3 --all2all yes"));

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 1 --num_vector_local 3 "
                    "--compute_method CPU --num_way 3",
                    "--num_proc_vector 1 --num_field 1 --num_vector_local 3 "
                    "--compute_method GPU --num_way 3 --all2all yes"));

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 1 --num_vector_local 18 "
                    "--compute_method CPU --num_way 3",
                    "--num_proc_vector 3 --num_field 1 --num_vector_local 6 "
                    "--compute_method CPU --num_way 3 --all2all yes"));

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 1 --num_vector_local 18 "
                    "--compute_method CPU --num_way 3",
                    "--num_proc_vector 3 --num_field 1 --num_vector_local 6 "
                    "--compute_method GPU --num_way 3 --all2all yes"));

  EXPECT_EQ(
      true,
      compare_3runs("--num_proc_vector 1 --num_field 2 --num_vector 16 "
                    "--compute_method REF --all2all yes --num_way 3",
                    "--num_proc_vector 3 --num_field 2 --num_vector 16 "
                    "--compute_method CPU --all2all yes --num_way 3",
                    "--num_proc_vector 3 --num_field 2 --num_vector 16 "
                    "--compute_method GPU --all2all yes --num_way 3"));

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_proc_field 1 "
                    "--num_field 13 --num_vector 6 "
                    "--compute_method REF --all2all yes --num_way 3",
                    "--num_proc_vector 1 --num_proc_field 5 "
                    "--num_field 13 --num_vector 6 "
                    "--compute_method GPU --all2all yes --num_way 3"));

  //----------
  //---3-way, all2all yes, large
  //----------

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method CPU --num_way 3",
                    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method CPU --num_way 3 --all2all yes"));

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method GPU --num_way 3",
                    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method GPU --num_way 3 --all2all yes"));

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method CPU --num_way 3 --all2all yes",
                    "--num_proc_vector 4 --num_field 100 --num_vector_local 12 "
                    "--compute_method CPU --num_way 3 --all2all yes"));

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method GPU --num_way 3 --all2all yes",
                    "--num_proc_vector 4 --num_field 100 --num_vector_local 12 "
                    "--compute_method GPU --num_way 3 --all2all yes"));

  //----------
  //---num_proc_field
  //----------

  EXPECT_EQ(true, compare_2runs("--num_proc_vector 1 --num_proc_field "
                                "1 --num_field 2 --num_vector_local 3 "
                                "--compute_method CPU --num_way 3",
                                "--num_proc_vector 1 --num_proc_field "
                                "2 --num_field 2 --num_vector_local 3 "
                                "--compute_method GPU --num_way 3"));

  EXPECT_EQ(true,
            compare_2runs("--num_proc_vector 1 --num_proc_field 1 --num_field "
                          "2 --num_vector_local 18"
                          " --compute_method CPU --num_way 3",
                          "--num_proc_vector 3 --num_proc_field 2 --num_field "
                          "2 --num_vector_local 6 "
                          " --compute_method GPU --num_way 3 --all2all yes"));
#endif

  //----------
  //---num_proc_repl, num_stage, 3-way
  //----------

  char options1[1024];
  char options2[1024];

  char options_template_1[] =
      "--metric_type czekanowski "
      "--num_field 4 --num_vector_local %i --compute_method %s --all2all yes "
      "--num_proc_vector %i --num_proc_repl %i "
      "--num_proc_field %i --num_way %i --num_stage %i";

  for (int gpu=0; gpu<=1; ++gpu) {
    for (int num_vector_local=6; num_vector_local<=18; num_vector_local+=12) {
      for (int num_proc_vector=1; num_proc_vector<=6; ++num_proc_vector) {
        for (int num_proc_repl=1; num_proc_repl<=6; ++num_proc_repl) {
          for (int num_stage=1; num_stage<=6; num_stage+=4) {
          const int num_proc_field = gpu ? 2 : 1;
            if (num_proc_vector * num_proc_field * num_proc_repl > PROCS_MAX) {
              continue;
            }
            const int num_way = 3;
            sprintf(options1, options_template_1,
                    num_vector_local*num_proc_vector,
                    "GPU", 1, 1,
                   1, num_way, 1);
            sprintf(options2, options_template_1, num_vector_local,
                    gpu ? "GPU" : "CPU", num_proc_vector, num_proc_repl,
                    num_proc_field, num_way, num_stage);
            test_2runs(options1, options2);
          }
        }
      }
    }
  }

  //----------
  //---num_proc_repl, num_phase, 3-way
  //----------

  char options_template_2[] =
      "--metric_type czekanowski "
      "--num_field 4 --num_vector_local %i --compute_method %s --all2all yes "
      "--num_proc_vector %i --num_proc_repl %i "
      "--num_proc_field %i --num_way %i --num_phase %i";

  for (int gpu=0; gpu<=1; ++gpu) {
    for (int num_vector_local=6; num_vector_local<=6; num_vector_local+=12) {
      for (int num_proc_vector=1; num_proc_vector<=4; ++num_proc_vector) {
        for (int num_proc_repl=1; num_proc_repl<=4; ++num_proc_repl) {
          const int npv = num_proc_vector;
          const int num_phase_max = 1==num_proc_repl ? npv*npv - 2*npv + 2 :
                                                      (npv+1)*(npv+2);
          const int num_phase_min = num_phase_max / 2;
          for (int num_phase=num_phase_min; num_phase<=num_phase_max;
               num_phase+=(num_phase_max-num_phase_min)) {
            if (num_phase < 1) {
              continue;
            }
            const int num_proc_field = gpu ? 2 : 1;
            if (num_proc_vector * num_proc_field * num_proc_repl > PROCS_MAX) {
              continue;
            }
            const int num_way = 3;
            sprintf(options1, options_template_2,
                    num_vector_local*num_proc_vector,
                    "GPU", 1, 1,
                   1, num_way, 1);
            sprintf(options2, options_template_2, num_vector_local,
                    gpu ? "GPU" : "CPU", num_proc_vector, num_proc_repl,
                    num_proc_field, num_way, num_phase);
            test_2runs(options1, options2);
          }
        }
      }
    }
  }

//  //----------
//  //---file output, 3-way
//  //----------
//
//  EXPECT_EQ(
//      true,
//      compare_2runs("--num_proc_field 1 --num_proc_vector 1 "
//                    "--num_field 7 --num_vector 18 "
//                    "--num_way 3 --metric_type czekanowski "
//                    "--compute_method REF --all2all yes",
//                    "--num_proc_field 1 --num_proc_vector 3 "
//                    "--num_field 7 --num_vector 18 "
//                    "--num_way 3 --metric_type czekanowski "
//                    "--compute_method GPU --all2all yes "
//                    "--verbosity 1 "
//                    "--output_file_stub test_czek_3way"));

} // DriverTest_czek3_

//=============================================================================

void DriverTest_ccc2_simple_compute_method(int compute_method) {

  using namespace comet;

  const int num_field = 5;
  const int num_vector = 2;

  std::string options = " --num_way " + std::to_string(2);
  options += " --metric_type " + std::string(MetricType::str(MetricType::CCC));
  options += " --all2all no";
  options += " --compute_method " + std::string(ComputeMethod::str(compute_method));
  options += " --num_proc_vector " + std::to_string(1);
  options += " --tc 4";
  options += " --verbosity 1";
  if (!can_run(options.c_str()))
    return;
  CEnv env_value(MPI_COMM_WORLD, options.c_str());
  CEnv* env = &env_value;

  GMDecompMgr dm(env_value);
  dm.allocate(true, false, num_field, num_vector, env->data_type_vectors());

  Vectors vectors(dm, env_value);
  vectors.allocate();
  vectors.initialize();

  if (env->is_proc_active()) {
    {
      const int G = 0;
      const int T = 1;
      int f = 0;
      const int i = 0;
      vectors.bits2_set(f++, i, 2 * G + 1 * T, *env);
      vectors.bits2_set(f++, i, 2 * T + 1 * T, *env);
      vectors.bits2_set(f++, i, 2 * T + 1 * T, *env);
      vectors.bits2_set(f++, i, 2 * T + 1 * T, *env);
      vectors.bits2_set(f++, i, 2 * T + 1 * T, *env);
    }
    {
      const int G = 0;
      const int A = 1;
      int f = 0;
      const int i = 1;
      vectors.bits2_set(f++, i, 2 * G + 1 * G, *env);
      vectors.bits2_set(f++, i, 2 * A + 1 * G, *env);
      vectors.bits2_set(f++, i, 2 * G + 1 * G, *env);
      vectors.bits2_set(f++, i, 2 * G + 1 * G, *env);
      vectors.bits2_set(f++, i, 2 * G + 1 * A, *env);
    }
  }

  GMMetrics metrics_value = GMMetrics_null();
  GMMetrics* metrics = &metrics_value;
  MetricsMem metrics_mem(env_value);
  GMMetrics_create(metrics, env->data_type_metrics(), &dm,
                   &metrics_mem, env);

  if (env->is_proc_active())
    printf("%s\n", options.c_str());

  ComputeMetrics::compute(*metrics, vectors, *env);

  Checksum cksum;
  Checksum cksum_local;
  Checksum::compute(cksum, cksum_local, *metrics, *env);
  Driver::print_output_sync(cksum, *env);

  if (env->is_proc_active()) {
    const double result00 = Metrics_ccc_duo_get_2(*metrics, 0, 0, 0, *env);
    const double result01 = Metrics_ccc_duo_get_2(*metrics, 0, 0, 1, *env);
    const double result10 = Metrics_ccc_duo_get_2(*metrics, 0, 1, 0, *env);
    const double result11 = Metrics_ccc_duo_get_2(*metrics, 0, 1, 1, *env);

    printf("COMPUTED: G G  %.5f\n", result00);
    printf("COMPUTED: G A  %.5f\n", result01);
    printf("COMPUTED: T G  %.5f\n", result10);
    printf("COMPUTED: T A  %.5f\n", result11);
    printf("\n");

    const double ref00 = .196;
    const double ref01 = .000;
    const double ref10 = .588;
    const double ref11 = .312;

    printf("EXPECTED: G G  %.5f\n", ref00);
    printf("EXPECTED: G A  %.5f\n", ref01);
    printf("EXPECTED: T G  %.5f\n", ref10);
    printf("EXPECTED: T A  %.5f\n", ref11);
    printf("\n");

    const double eps = 1.e-5;

    EXPECT_EQ(true, fabs(result00 - ref00) < eps);
    EXPECT_EQ(true, fabs(result01 - ref01) < eps);
    EXPECT_EQ(true, fabs(result10 - ref10) < eps);
    EXPECT_EQ(true, fabs(result11 - ref11) < eps);
  }

  GMMetrics_destroy(metrics, env);
  vectors.deallocate();
  dm.deallocate();
} // DriverTest_ccc2_simple_compute_method

//=============================================================================

void DriverTest_ccc2_simple_() {
  DriverTest_ccc2_simple_compute_method(comet::ComputeMethod::REF);
  DriverTest_ccc2_simple_compute_method(comet::ComputeMethod::CPU);
  DriverTest_ccc2_simple_compute_method(comet::ComputeMethod::GPU);
}

//=============================================================================

void DriverTest_ccc2_simple_sparse_compute_method(int compute_method) {

  using namespace comet;

  const int num_field = 5;
  const int num_vector = 2;

  std::string options = " --num_way " + std::to_string(2);
  options += " --metric_type " + std::string(MetricType::str(MetricType::CCC));
  options += " --all2all no";
  options += " --sparse yes";
  options += " --compute_method " + std::string(ComputeMethod::str(compute_method));
  options += " --num_proc_vector " + std::to_string(1);
  options += " --tc 4";
  options += " --verbosity 1";
  if (!can_run(options.c_str()))
    return;
  CEnv env_value(MPI_COMM_WORLD, options.c_str());
  CEnv* env = &env_value;

  GMDecompMgr dm(env_value);
  dm.allocate(true, false, num_field, num_vector, env->data_type_vectors());

  Vectors vectors(dm, env_value);
  vectors.allocate();
  vectors.initialize();

  if (env->is_proc_active()) {
    const int UN = 2 * 1 + 1 * 0;
    {
      const int G = 0;
      const int T = 1;
    //const int GG =  2 * G + 1 * G
      const int GT =  2 * G + 1 * T;
    //const int TG =  2 * G + 1 * T
      const int TT =  2 * T + 1 * T;
      int f = 0;
      const int i = 0;
      vectors.bits2_set(f++, i, GT, *env);
      vectors.bits2_set(f++, i, TT, *env);
      vectors.bits2_set(f++, i, TT, *env);
      vectors.bits2_set(f++, i, UN, *env);
      vectors.bits2_set(f++, i, UN, *env);
    }
    {
      const int G = 0;
      const int A = 1;
      const int GG =  2 * G + 1 * G;
      const int GA =  2 * G + 1 * A;
      const int AG =  2 * G + 1 * A;
    //const int AA =  2 * Q + 1 * A
      int f = 0;
      const int i = 1;
      vectors.bits2_set(f++, i, GG, *env);
      vectors.bits2_set(f++, i, AG, *env);
      vectors.bits2_set(f++, i, GG, *env);
      vectors.bits2_set(f++, i, UN, *env);
      vectors.bits2_set(f++, i, GA, *env);
    }
  }

  GMMetrics metrics_value = GMMetrics_null();
  GMMetrics* metrics = &metrics_value;
  MetricsMem metrics_mem(env_value);
  GMMetrics_create(metrics, env->data_type_metrics(), &dm,
                   &metrics_mem, env);

  if (env->is_proc_active())
    printf("%s\n", options.c_str());

  ComputeMetrics::compute(*metrics, vectors, *env);

  Checksum cksum;
  Checksum cksum_local;
  Checksum::compute(cksum, cksum_local, *metrics, *env);
  Driver::print_output_sync(cksum, *env);

  if (env->is_proc_active()) {
    const double result00 = Metrics_ccc_duo_get_2(*metrics, 0, 0, 0, *env);
    const double result01 = Metrics_ccc_duo_get_2(*metrics, 0, 0, 1, *env);
    const double result10 = Metrics_ccc_duo_get_2(*metrics, 0, 1, 0, *env);
    const double result11 = Metrics_ccc_duo_get_2(*metrics, 0, 1, 1, *env);

    printf("COMPUTED: G G  %.5f\n", result00);
    printf("COMPUTED: G A  %.5f\n", result01);
    printf("COMPUTED: T G  %.5f\n", result10);
    printf("COMPUTED: T A  %.5f\n", result11);
    printf("\n");

    const double s0_0 = 1;
    const double s0_1 = 5;

    const double s1_0 = 6;
    const double s1_1 = 2;

    const double c0 = s0_0 + s0_1;
    const double c1 = s1_0 + s1_1;

    const double f0_0 = s0_0 / c0;
    const double f0_1 = s0_1 / c0;

    const double f1_0 = s1_0 / c1;
    const double f1_1 = s1_1 / c1;

    const double r0_00 = 2;
    const double r0_01 = 0;
    const double r0_10 = 2;
    const double r0_11 = 0;

    const double r1_00 = 0;
    const double r1_01 = 0;
    const double r1_10 = 2;
    const double r1_11 = 2;

    const double r2_00 = 0;
    const double r2_01 = 0;
    const double r2_10 = 4;
    const double r2_11 = 0;

    //const double r3_*** = 0
    //const double r4_*** = 0

    const double r_00 = r0_00 + r1_00 + r2_00;
    const double r_01 = r0_01 + r1_01 + r2_01;
    const double r_10 = r0_10 + r1_10 + r2_10;
    const double r_11 = r0_11 + r1_11 + r2_11;

    const double c = r_00 + r_01 + r_10 + r_11;

    const double f_00 = r_00 / c;
    const double f_01 = r_01 / c;
    const double f_10 = r_10 / c;
    const double f_11 = r_11 / c;

    const double fm = 9 / (double) 2;
    const double cp = 2 / (double) 3;

    const double ref00 = fm * f_00 * ( 1 - cp * f0_0 ) * ( 1 - cp * f1_0 );
    const double ref01 = fm * f_01 * ( 1 - cp * f0_0 ) * ( 1 - cp * f1_1 );
    const double ref10 = fm * f_10 * ( 1 - cp * f0_1 ) * ( 1 - cp * f1_0 );
    const double ref11 = fm * f_11 * ( 1 - cp * f0_1 ) * ( 1 - cp * f1_1 );

    printf("EXPECTED: G G  %.5f\n", ref00);
    printf("EXPECTED: G A  %.5f\n", ref01);
    printf("EXPECTED: T G  %.5f\n", ref10);
    printf("EXPECTED: T A  %.5f\n", ref11);
    printf("\n");

    const double eps = 1.e-5;

    EXPECT_EQ(true, fabs(result00 - ref00) < eps);
    EXPECT_EQ(true, fabs(result01 - ref01) < eps);
    EXPECT_EQ(true, fabs(result10 - ref10) < eps);
    EXPECT_EQ(true, fabs(result11 - ref11) < eps);
  }

  GMMetrics_destroy(metrics, env);
  vectors.deallocate();
  dm.deallocate();
} // DriverTest_ccc2_simple_sparse_compute_method

//=============================================================================

void DriverTest_ccc2_simple_sparse_() {
  DriverTest_ccc2_simple_sparse_compute_method(comet::ComputeMethod::REF);
  DriverTest_ccc2_simple_sparse_compute_method(comet::ComputeMethod::CPU);
  DriverTest_ccc2_simple_sparse_compute_method(comet::ComputeMethod::GPU);
}

//=============================================================================

void DriverTest_duo2_simple_sparse_(int compute_method, double tLL, double tLH,
                                    double tHH, double tLLHH) {

  using namespace comet;

  const int num_field = 5;
  const int num_vector = 2;

  std::string options = " --num_way " + std::to_string(2);
  options += " --metric_type " + std::string(MetricType::str(MetricType::DUO));
  options += " --all2all no";
  options += " --sparse yes";
  options += " --compute_method " + std::string(ComputeMethod::str(compute_method));
  options += " --num_proc_vector " + std::to_string(1);
  options += " --tc 4";
  options += " --verbosity 1";
  options += " --threshold " + std::to_string(tLL) + "," + std::to_string(tLH) + ","
                             + std::to_string(tHH) + "," + std::to_string(tLLHH);
  if (!can_run(options.c_str()))
    return;
  CEnv env_value(MPI_COMM_WORLD, options.c_str());
  CEnv* env = &env_value;

  GMDecompMgr dm(env_value);
  dm.allocate(true, false, num_field, num_vector, env->data_type_vectors());

  Vectors vectors(dm, env_value);
  vectors.allocate();
  vectors.initialize();

  if (env->is_proc_active()) {
    // vector entry choices, binary representation
    const int MIN = 2 * (0) + 1 * (0);
    const int MAX = 2 * (1) + 1 * (1);
    const int UNK = 2 * (1) + 1 * (0);
    // define first vector
    {
      int f = 0;
      const int i = 0;
      vectors.bits2_set(f++, i, MIN, *env);
      vectors.bits2_set(f++, i, MAX, *env);
      vectors.bits2_set(f++, i, MIN, *env);
      vectors.bits2_set(f++, i, UNK, *env);
      vectors.bits2_set(f++, i, UNK, *env);
    }
    // define second vector
    {
      int f = 0;
      const int i = 1;
      vectors.bits2_set(f++, i, UNK, *env);
      vectors.bits2_set(f++, i, MAX, *env);
      vectors.bits2_set(f++, i, MAX, *env);
      vectors.bits2_set(f++, i, MAX, *env);
      vectors.bits2_set(f++, i, UNK, *env);
    }
  }

  GMMetrics metrics_value = GMMetrics_null();
  GMMetrics* metrics = &metrics_value;
  MetricsMem metrics_mem(env_value);
  GMMetrics_create(metrics, env->data_type_metrics(), &dm,
                   &metrics_mem, env);

  if (env->is_proc_active())
    printf("%s\n", options.c_str());

  ComputeMetrics::compute(*metrics, vectors, *env);

  Checksum cksum;
  Checksum cksum_local;
  Checksum::compute(cksum, cksum_local, *metrics, *env);
  Driver::print_output_sync(cksum, *env);

  if (env->is_proc_active()) {

    const double result00 =
      Metrics_is_pass_threshold_noshrink(*metrics, 0, 0, 0, *env) ?
      Metrics_ccc_duo_get_2(*metrics, 0, 0, 0, *env) : 0;
    const double result01 =
      Metrics_is_pass_threshold_noshrink(*metrics, 0, 0, 1, *env) ?
      Metrics_ccc_duo_get_2(*metrics, 0, 0, 1, *env) : 0;
    const double result10 =
      Metrics_is_pass_threshold_noshrink(*metrics, 0, 1, 0, *env) ?
      Metrics_ccc_duo_get_2(*metrics, 0, 1, 0, *env) : 0;
    const double result11 =
      Metrics_is_pass_threshold_noshrink(*metrics, 0, 1, 1, *env) ?
      Metrics_ccc_duo_get_2(*metrics, 0, 1, 1, *env) : 0;

    //const double result00 = Metrics_ccc_duo_get_2(*metrics, 0, 0, 0, *env);
    //const double result01 = Metrics_ccc_duo_get_2(*metrics, 0, 0, 1, *env);
    //const double result10 = Metrics_ccc_duo_get_2(*metrics, 0, 1, 0, *env);
    //const double result11 = Metrics_ccc_duo_get_2(*metrics, 0, 1, 1, *env);

    printf("COMPUTED: MIN MIN  %.17f\n", result00);
    printf("COMPUTED: MIN MAX  %.17f\n", result01);
    printf("COMPUTED: MAX MIN  %.17f\n", result10);
    printf("COMPUTED: MAX MAX  %.17f\n", result11);
    printf("\n");

    // calculate by hand the expected DUO value.

    // recall that num_field = 5

    const double s0_0 = 2; // vector 0, number of MINs
    const double s0_1 = 1; // vector 0, number of MAXs

    const double s1_0 = 0;
    const double s1_1 = 3;

    const double c0 = s0_0 + s0_1; // vector 0, number of MINs and MAXs
    const double c1 = s1_0 + s1_1;

    // const double unk_0 = num_field - c0; // = 2 // for vector 0, num UNK
    // const double unk_1 = num_field - c1; // = 2 // for vector 1, num UNK

    const double f0_0 = s0_0 / c0; // vector 0, f_0(MIN)
    const double f0_1 = s0_1 / c0; // vector 0, f_0(MAX)

    const double f1_0 = s1_0 / c1;
    const double f1_1 = s1_1 / c1;

    // numerators from table computation
    const double r_00 = 0;
    const double r_01 = 1;
    const double r_10 = 0;
    const double r_11 = 1;

    // sum of all table entries
    const double c = r_00 + r_01 + r_10 + r_11;

    // Dij of DUO
    const double d_00 = r_00 / c;
    const double d_01 = r_01 / c;
    const double d_10 = r_10 / c;
    const double d_11 = r_11 / c;

    // Constants needed by method
    const double fm = 4;
    const double cp = 2 / (double) 3;

    // DUO values
    const double ref00 = fm * d_00 * ( 1 - cp * f0_0 ) * ( 1 - cp * f1_0 );
    const double ref01 = fm * d_01 * ( 1 - cp * f0_0 ) * ( 1 - cp * f1_1 );
    const double ref10 = fm * d_10 * ( 1 - cp * f0_1 ) * ( 1 - cp * f1_0 );
    const double ref11 = fm * d_11 * ( 1 - cp * f0_1 ) * ( 1 - cp * f1_1 );

    const double ref00t = ref00 > tLL || (ref00 + ref11 > tLLHH &&
      ref00 > 0 && ref11 > 0) ? ref00 : (double)0;
    const double ref01t = ref01 > tLH                  ? ref01 : (double)0;
    const double ref10t = ref01 > tLH                  ? ref10 : (double)0;
    const double ref11t = ref11 > tHH || (ref00 + ref11 > tLLHH &&
      ref00 > 0 && ref11 > 0) ? ref11 : (double)0;

    printf("EXPECTED: MIN MIN  %.17f\n", ref00t);
    printf("EXPECTED: MIN MAX  %.17f\n", ref01t);
    printf("EXPECTED: MAX MIN  %.17f\n", ref10t);
    printf("EXPECTED: MAX MAX  %.17f\n", ref11t);
    printf("\n");

    const double eps = 1.e-5;

    EXPECT_EQ(true, fabs(result00 - ref00t) < eps);
    EXPECT_EQ(true, fabs(result01 - ref01t) < eps);
    EXPECT_EQ(true, fabs(result10 - ref10t) < eps);
    EXPECT_EQ(true, fabs(result11 - ref11t) < eps);
  }

  GMMetrics_destroy(metrics, env);
  vectors.deallocate();
  dm.deallocate();
} // DriverTest_duo2_simple_sparse_

//=============================================================================

void DriverTest_duo2_simple_sparse_(int compute_method) {

//  if (comet::ComputeMethod::REF == compute_method ||
//      comet::BuildHas::DOUBLE_PREC) {
//
//    // REF method and DP build do not remove threshold-failed values,
//    // so don't test that.
//    const double tLL = .01;
//    const double tLH = .01;
//    const double tHH = .01;
//    const double tLLHH = .01;
//
//    DriverTest_duo2_simple_sparse_(compute_method, tLL, tLH, tHH, tLLHH);
//
//  } else {

    const double tLL = .01;
    for (double tLH : {.37, .38}) {
    for (double tHH : {.51, .52}) {
    for (double tLLHH : {.51, .52}) {
      DriverTest_duo2_simple_sparse_(compute_method, tLL, tLH, tHH, tLLHH);
    }
    }
    }

// }
}

//=============================================================================

void DriverTest_duo2_simple_sparse_() {
  DriverTest_duo2_simple_sparse_(comet::ComputeMethod::REF);
  DriverTest_duo2_simple_sparse_(comet::ComputeMethod::CPU);
  DriverTest_duo2_simple_sparse_(comet::ComputeMethod::GPU);
}

//=============================================================================

void DriverTest_ccc3_simple_compute_method(int compute_method) {

  using namespace comet;

  const int num_field = 10;
  const int num_vector = 3;

  std::string options = " --num_way " + std::to_string(3);
  options += " --metric_type " + std::string(MetricType::str(MetricType::CCC));
  options += " --all2all yes";
  options += " --compute_method " + std::string(ComputeMethod::str(compute_method));
  options += " --num_proc_vector " + std::to_string(1);
  options += " --tc 4";
  options += " --verbosity 1";
  if (!can_run(options.c_str()))
    return;
  CEnv env_value(MPI_COMM_WORLD, options.c_str());
  CEnv* env = &env_value;

  GMDecompMgr dm(env_value);
  dm.allocate(true, false, num_field, num_vector, env->data_type_vectors());

  Vectors vectors(dm, env_value);
  vectors.allocate();
  vectors.initialize();

  if (env->is_proc_active()) {
    {
      const int A = 0;
      const int T = 1;
      int f = 0;
      const int i = 0;
      vectors.bits2_set(f++, i, 2 * A + 1 * A, *env);
      vectors.bits2_set(f++, i, 2 * A + 1 * T, *env);
      vectors.bits2_set(f++, i, 2 * T + 1 * T, *env);
      vectors.bits2_set(f++, i, 2 * A + 1 * T, *env);
      vectors.bits2_set(f++, i, 2 * A + 1 * T, *env);
      vectors.bits2_set(f++, i, 2 * A + 1 * T, *env);
      vectors.bits2_set(f++, i, 2 * A + 1 * A, *env);
      vectors.bits2_set(f++, i, 2 * A + 1 * A, *env);
      vectors.bits2_set(f++, i, 2 * A + 1 * A, *env);
      vectors.bits2_set(f++, i, 2 * A + 1 * A, *env);
    }
    {
      const int A = 0;
      const int T = 1;
      int f = 0;
      const int i = 1;
      vectors.bits2_set(f++, i, 2 * A + 1 * T, *env);
      vectors.bits2_set(f++, i, 2 * A + 1 * A, *env);
      vectors.bits2_set(f++, i, 2 * A + 1 * A, *env);
      vectors.bits2_set(f++, i, 2 * A + 1 * T, *env);
      vectors.bits2_set(f++, i, 2 * A + 1 * T, *env);
      vectors.bits2_set(f++, i, 2 * A + 1 * A, *env);
      vectors.bits2_set(f++, i, 2 * A + 1 * T, *env);
      vectors.bits2_set(f++, i, 2 * A + 1 * A, *env);
      vectors.bits2_set(f++, i, 2 * A + 1 * A, *env);
      vectors.bits2_set(f++, i, 2 * A + 1 * A, *env);
    }
    {
      const int A = 0;
      const int T = 1;
      int f = 0;
      const int i = 2;
      vectors.bits2_set(f++, i, 2 * A + 1 * A, *env);
      vectors.bits2_set(f++, i, 2 * A + 1 * A, *env);
      vectors.bits2_set(f++, i, 2 * A + 1 * A, *env);
      vectors.bits2_set(f++, i, 2 * A + 1 * T, *env);
      vectors.bits2_set(f++, i, 2 * A + 1 * T, *env);
      vectors.bits2_set(f++, i, 2 * A + 1 * A, *env);
      vectors.bits2_set(f++, i, 2 * T + 1 * T, *env);
      vectors.bits2_set(f++, i, 2 * A + 1 * A, *env);
      vectors.bits2_set(f++, i, 2 * A + 1 * A, *env);
      vectors.bits2_set(f++, i, 2 * A + 1 * A, *env);
    }
  }

  GMMetrics metrics_value = GMMetrics_null();
  GMMetrics* metrics = &metrics_value;
  MetricsMem metrics_mem(env_value);
  GMMetrics_create(metrics, env->data_type_metrics(), &dm,
                   &metrics_mem, env);

  if (env->is_proc_active())
    printf("%s\n", options.c_str());

  ComputeMetrics::compute(*metrics, vectors, *env);

  Checksum cksum;
  Checksum cksum_local;
  Checksum::compute(cksum, cksum_local, *metrics, *env);
  Driver::print_output_sync(cksum, *env);

  if (env->is_proc_active()) {
    const double result000 = Metrics_ccc_duo_get_3(*metrics, 0, 0, 0, 0, *env);
    const double result001 = Metrics_ccc_duo_get_3(*metrics, 0, 0, 0, 1, *env);
    const double result010 = Metrics_ccc_duo_get_3(*metrics, 0, 0, 1, 0, *env);
    const double result011 = Metrics_ccc_duo_get_3(*metrics, 0, 0, 1, 1, *env);
    const double result100 = Metrics_ccc_duo_get_3(*metrics, 0, 1, 0, 0, *env);
    const double result101 = Metrics_ccc_duo_get_3(*metrics, 0, 1, 0, 1, *env);
    const double result110 = Metrics_ccc_duo_get_3(*metrics, 0, 1, 1, 0, *env);
    const double result111 = Metrics_ccc_duo_get_3(*metrics, 0, 1, 1, 1, *env);

    printf("COMPUTED: A A A  %.8f\n", result000);
    printf("COMPUTED: A A T  %.5f\n", result001);
    printf("COMPUTED: A T A  %.8f\n", result010);
    printf("COMPUTED: A T T  %.8f\n", result011);
    printf("COMPUTED: T A A  %.8f\n", result100);
    printf("COMPUTED: T A T  %.8f\n", result101);
    printf("COMPUTED: T T A  %.8f\n", result110);
    printf("COMPUTED: T T T  %.8f\n", result111);
    printf("\n");

    const double fm = 9 / (double) 2;
    const double ref000 = fm * .055;
    const double ref001 = fm * .016;
    const double ref010 = fm * .016;
    const double ref011 = fm * .030;
    const double ref100 = fm * .039;
    const double ref101 = fm * .008;
    const double ref110 = fm * .008;
    const double ref111 = fm * .015;

    printf("EXPECTED: A A A  %.8f\n", ref000);
    printf("EXPECTED: A A T  %.5f\n", ref001);
    printf("EXPECTED: A T A  %.8f\n", ref010);
    printf("EXPECTED: A T T  %.8f\n", ref011);
    printf("EXPECTED: T A A  %.8f\n", ref100);
    printf("EXPECTED: T A T  %.8f\n", ref101);
    printf("EXPECTED: T T A  %.8f\n", ref110);
    printf("EXPECTED: T T T  %.8f\n", ref111);
    printf("\n");

    const double eps = 1.e-3;

    EXPECT_EQ(true, fabs(result000 - ref000) < eps);
    EXPECT_EQ(true, fabs(result001 - ref001) < eps);
    EXPECT_EQ(true, fabs(result010 - ref010) < eps);
    EXPECT_EQ(true, fabs(result011 - ref011) < eps);
    EXPECT_EQ(true, fabs(result100 - ref100) < eps);
    EXPECT_EQ(true, fabs(result101 - ref101) < eps);
    EXPECT_EQ(true, fabs(result110 - ref110) < eps);
    EXPECT_EQ(true, fabs(result111 - ref111) < eps);
  }

  GMMetrics_destroy(metrics, env);
  vectors.deallocate();
  dm.deallocate();
} // DriverTest_ccc3_simple_compute_method

//=============================================================================

void DriverTest_ccc3_simple_() {
  DriverTest_ccc3_simple_compute_method(comet::ComputeMethod::REF);
  DriverTest_ccc3_simple_compute_method(comet::ComputeMethod::CPU);
  DriverTest_ccc3_simple_compute_method(comet::ComputeMethod::GPU);
}

//=============================================================================

void DriverTest_ccc3_simple_sparse_compute_method(int compute_method) {

  using namespace comet;

  const int num_field = 10;
  const int num_vector = 3;

  std::string options = " --num_way " + std::to_string(3);
  options += " --metric_type " + std::string(MetricType::str(MetricType::CCC));
  options += " --all2all yes";
  options += " --sparse yes";
  options += " --compute_method " + std::string(ComputeMethod::str(compute_method));
  options += " --num_proc_vector " + std::to_string(1);
  options += " --tc 4";
  options += " --verbosity 1";
  options += " --threshold .0025";
  if (!can_run(options.c_str()))
    return;
  CEnv env_value(MPI_COMM_WORLD, options.c_str());
  CEnv* env = &env_value;

  GMDecompMgr dm(env_value);
  dm.allocate(true, false, num_field, num_vector, env->data_type_vectors());

  Vectors vectors(dm, env_value);
  vectors.allocate();
  vectors.initialize();

  if (env->is_proc_active()) {
    const int UN = 2 * 1 + 1 * 0;
    {
      const int A = 0;
      const int T = 1;
      const int AA =  2 * A + 1 * A;
      const int AT =  2 * A + 1 * T;
    //const int TA =  2 * A + 1 * T
      const int TT =  2 * T + 1 * T;
      int f = 0;
      const int i = 0;
      vectors.bits2_set(f++, i, AA, *env);
      vectors.bits2_set(f++, i, AT, *env);
      vectors.bits2_set(f++, i, TT, *env);
      vectors.bits2_set(f++, i, UN, *env);
      vectors.bits2_set(f++, i, AT, *env);
      vectors.bits2_set(f++, i, AT, *env);
      vectors.bits2_set(f++, i, UN, *env);
      vectors.bits2_set(f++, i, UN, *env);
      vectors.bits2_set(f++, i, AA, *env);
      vectors.bits2_set(f++, i, AA, *env);
    }
    {
      const int A = 0;
      const int T = 1;
      const int AA =  2 * A + 1 * A;
      const int AT =  2 * A + 1 * T;
    //const int TA =  2 * A + 1 * T
    //const int TT =  2 * T + 1 * T
      int f = 0;
      const int i = 1;
      vectors.bits2_set(f++, i, AT, *env);
      vectors.bits2_set(f++, i, AA, *env);
      vectors.bits2_set(f++, i, AA, *env);
      vectors.bits2_set(f++, i, AT, *env);
      vectors.bits2_set(f++, i, AT, *env);
      vectors.bits2_set(f++, i, AA, *env);
      vectors.bits2_set(f++, i, AT, *env);
      vectors.bits2_set(f++, i, UN, *env);
      vectors.bits2_set(f++, i, UN, *env);
      vectors.bits2_set(f++, i, AA, *env);
    }
    {
      const int A = 0;
      const int T = 1;
      const int AA =  2 * A + 1 * A;
      const int AT =  2 * A + 1 * T;
    //const int TA =  2 * A + 1 * T
    //const int TT =  2 * T + 1 * T
      int f = 0;
      const int i = 2;
      vectors.bits2_set(f++, i, AA, *env);
      vectors.bits2_set(f++, i, AA, *env);
      vectors.bits2_set(f++, i, AA, *env);
      vectors.bits2_set(f++, i, AT, *env);
      vectors.bits2_set(f++, i, AT, *env);
      vectors.bits2_set(f++, i, AA, *env);
      vectors.bits2_set(f++, i, UN, *env);
      vectors.bits2_set(f++, i, AA, *env);
      vectors.bits2_set(f++, i, AA, *env);
      vectors.bits2_set(f++, i, UN, *env);
    }
  }

  GMMetrics metrics_value = GMMetrics_null();
  GMMetrics* metrics = &metrics_value;
  MetricsMem metrics_mem(env_value);
  GMMetrics_create(metrics, env->data_type_metrics(), &dm,
                   &metrics_mem, env);

  if (env->is_proc_active())
    printf("%s\n", options.c_str());

  ComputeMetrics::compute(*metrics, vectors, *env);

  Checksum cksum;
  Checksum cksum_local;
  Checksum::compute(cksum, cksum_local, *metrics, *env);
  Driver::print_output_sync(cksum, *env);

  if (env->is_proc_active()) {
    const double result000 = Metrics_ccc_duo_get_3(*metrics, 0, 0, 0, 0, *env);
    const double result001 = Metrics_ccc_duo_get_3(*metrics, 0, 0, 0, 1, *env);
    const double result010 = Metrics_ccc_duo_get_3(*metrics, 0, 0, 1, 0, *env);
    const double result011 = Metrics_ccc_duo_get_3(*metrics, 0, 0, 1, 1, *env);
    const double result100 = Metrics_ccc_duo_get_3(*metrics, 0, 1, 0, 0, *env);
    const double result101 = Metrics_ccc_duo_get_3(*metrics, 0, 1, 0, 1, *env);
    const double result110 = Metrics_ccc_duo_get_3(*metrics, 0, 1, 1, 0, *env);
    const double result111 = Metrics_ccc_duo_get_3(*metrics, 0, 1, 1, 1, *env);

    printf("COMPUTED: A A A  %.8f\n", result000);
    printf("COMPUTED: A A T  %.5f\n", result001);
    printf("COMPUTED: A T A  %.8f\n", result010);
    printf("COMPUTED: A T T  %.8f\n", result011);
    printf("COMPUTED: T A A  %.8f\n", result100);
    printf("COMPUTED: T A T  %.8f\n", result101);
    printf("COMPUTED: T T A  %.8f\n", result110);
    printf("COMPUTED: T T T  %.8f\n", result111);
    printf("\n");

// TODO: here and elsewhere, c0 -> cX2_ etc.
    const double s0_0 = 9;
    const double s0_1 = 5;
    const double c0 = s0_0 + s0_1;

    const double s1_0 = 12;
    const double s1_1 = 4;
    const double c1 = s1_0 + s1_1;

    const double s2_0 = 14;
    const double s2_1 = 2;
    const double c2 = s2_0 + s2_1;

    const double f0_0 = s0_0 / c0;
    const double f0_1 = s0_1 / c0;
    const double f1_0 = s1_0 / c1;
    const double f1_1 = s1_1 / c1;
    const double f2_0 = s2_0 / c2;
    const double f2_1 = s2_1 / c2;

    const double r0_000 = 4;
    const double r0_001 = 0;
    const double r0_010 = 4;
    const double r0_011 = 0;
    const double r0_100 = 0;
    const double r0_101 = 0;
    const double r0_110 = 0;
    const double r0_111 = 0;

    const double r1_000 = 4;
    const double r1_001 = 0;
    const double r1_010 = 0;
    const double r1_011 = 0;
    const double r1_100 = 4;
    const double r1_101 = 0;
    const double r1_110 = 0;
    const double r1_111 = 0;

    const double r2_000 = 0;
    const double r2_001 = 0;
    const double r2_010 = 0;
    const double r2_011 = 0;
    const double r2_100 = 8;
    const double r2_101 = 0;
    const double r2_110 = 0;
    const double r2_111 = 0;

    const double r4_000 = 1;
    const double r4_001 = 1;
    const double r4_010 = 1;
    const double r4_011 = 1;
    const double r4_100 = 1;
    const double r4_101 = 1;
    const double r4_110 = 1;
    const double r4_111 = 1;

    const double r5_000 = 4;
    const double r5_001 = 0;
    const double r5_010 = 0;
    const double r5_011 = 0;
    const double r5_100 = 4;
    const double r5_101 = 0;
    const double r5_110 = 0;
    const double r5_111 = 0;

    //const double r3_*** = 0
    //const double r6_*** = 0
    //const double r7_*** = 0
    //const double r8_*** = 0
    //const double r9_*** = 0

    const double r_000 = r0_000 + r1_000 + r2_000 + r4_000 + r5_000;
    const double r_001 = r0_001 + r1_001 + r2_001 + r4_001 + r5_001;
    const double r_010 = r0_010 + r1_010 + r2_010 + r4_010 + r5_010;
    const double r_011 = r0_011 + r1_011 + r2_011 + r4_011 + r5_011;
    const double r_100 = r0_100 + r1_100 + r2_100 + r4_100 + r5_100;
    const double r_101 = r0_101 + r1_101 + r2_101 + r4_101 + r5_101;
    const double r_110 = r0_110 + r1_110 + r2_110 + r4_110 + r5_110;
    const double r_111 = r0_111 + r1_111 + r2_111 + r4_111 + r5_111;

    const double c = r_000 + r_001 + r_010 + r_011 +
                     r_100 + r_101 + r_110 + r_111;

    const double f_000 = r_000 / c;
    const double f_001 = r_001 / c;
    const double f_010 = r_010 / c;
    const double f_011 = r_011 / c;
    const double f_100 = r_100 / c;
    const double f_101 = r_101 / c;
    const double f_110 = r_110 / c;
    const double f_111 = r_111 / c;

    const double fm = 9 / (double) 2;
    const double cp = 2 / (double) 3;

    const double ref000 = fm * f_000 * (1-cp*f0_0) * (1-cp*f1_0) * (1-cp*f2_0);
    const double ref001 = fm * f_001 * (1-cp*f0_0) * (1-cp*f1_0) * (1-cp*f2_1);
    const double ref010 = fm * f_010 * (1-cp*f0_0) * (1-cp*f1_1) * (1-cp*f2_0);
    const double ref011 = fm * f_011 * (1-cp*f0_0) * (1-cp*f1_1) * (1-cp*f2_1);
    const double ref100 = fm * f_100 * (1-cp*f0_1) * (1-cp*f1_0) * (1-cp*f2_0);
    const double ref101 = fm * f_101 * (1-cp*f0_1) * (1-cp*f1_0) * (1-cp*f2_1);
    const double ref110 = fm * f_110 * (1-cp*f0_1) * (1-cp*f1_1) * (1-cp*f2_0);
    const double ref111 = fm * f_111 * (1-cp*f0_1) * (1-cp*f1_1) * (1-cp*f2_1);

    printf("EXPECTED: A A A  %.8f\n", ref000);
    printf("EXPECTED: A A T  %.5f\n", ref001);
    printf("EXPECTED: A T A  %.8f\n", ref010);
    printf("EXPECTED: A T T  %.8f\n", ref011);
    printf("EXPECTED: T A A  %.8f\n", ref100);
    printf("EXPECTED: T A T  %.8f\n", ref101);
    printf("EXPECTED: T T A  %.8f\n", ref110);
    printf("EXPECTED: T T T  %.8f\n", ref111);
    printf("\n");

    //const double ref000 = .055
    //const double ref001 = .016
    //const double ref010 = .016
    //const double ref011 = .030
    //const double ref100 = .039
    //const double ref101 = .008
    //const double ref110 = .008
    //const double ref111 = .015

    const double eps = 1.e-3;

    EXPECT_EQ(true, fabs(result000 - ref000) < eps);
    EXPECT_EQ(true, fabs(result001 - ref001) < eps);
    EXPECT_EQ(true, fabs(result010 - ref010) < eps);
    EXPECT_EQ(true, fabs(result011 - ref011) < eps);
    EXPECT_EQ(true, fabs(result100 - ref100) < eps);
    EXPECT_EQ(true, fabs(result101 - ref101) < eps);
    EXPECT_EQ(true, fabs(result110 - ref110) < eps);
    EXPECT_EQ(true, fabs(result111 - ref111) < eps);
  }

  GMMetrics_destroy(metrics, env);
  vectors.deallocate();
  dm.deallocate();
} // DriverTest_ccc3_simple_sparse_compute_method

//=============================================================================

void DriverTest_ccc3_simple_sparse_() {
  DriverTest_ccc3_simple_sparse_compute_method(comet::ComputeMethod::REF);
  DriverTest_ccc3_simple_sparse_compute_method(comet::ComputeMethod::CPU);
  DriverTest_ccc3_simple_sparse_compute_method(comet::ComputeMethod::GPU);
}

//=============================================================================

void DriverTest_duo3_simple_sparse_(int compute_method, double tLLL,
  double tLLH, double tLHH, double tHHH, double tLLLHHH) {

  using namespace comet;

  const int num_field = 10;
  const int num_vector = 3;

  std::string options = " --num_way " + std::to_string(3);
  options += " --metric_type " + std::string(MetricType::str(MetricType::DUO));
  options += " --all2all yes";
  options += " --sparse yes";
  options += " --compute_method " + std::string(ComputeMethod::str(compute_method));
  options += " --num_proc_vector " + std::to_string(1);
  options += " --tc 4";
  options += " --verbosity 1";
  options += " --threshold " + std::to_string(tLLL) + "," + std::to_string(tLLH) + ","
    + std::to_string(tLHH) + "," + std::to_string(tHHH) + "," + std::to_string(tLLLHHH);
  if (!can_run(options.c_str()))
    return;
  CEnv env_value(MPI_COMM_WORLD, options.c_str());
  CEnv* env = &env_value;

  if (ComputeMethod::GPU == compute_method && env->tc_eff() == TC::NO) {
    // DUO 3-way GPU requires tc.
    return;
  }

  GMDecompMgr dm(env_value);
  dm.allocate(true, false, num_field, num_vector, env->data_type_vectors());

  Vectors vectors(dm, env_value);
  vectors.allocate();
  vectors.initialize();

  if (env->is_proc_active()) {
    const int UN = 2 * 1 + 1 * 0;
    {
      const int _A = 0;
      const int _T = 1;
      int f = 0;
      const int i = 0;
      vectors.bits2_set(f++, i, _A, *env);
      vectors.bits2_set(f++, i, _T, *env);
      vectors.bits2_set(f++, i, _T, *env);
      vectors.bits2_set(f++, i, UN, *env);
      vectors.bits2_set(f++, i, _T, *env);
      vectors.bits2_set(f++, i, _T, *env);
      vectors.bits2_set(f++, i, UN, *env);
      vectors.bits2_set(f++, i, UN, *env);
      vectors.bits2_set(f++, i, _A, *env);
      vectors.bits2_set(f++, i, _A, *env);
    }
    {
      const int _A = 0;
      const int _T = 1;
      int f = 0;
      const int i = 1;
      vectors.bits2_set(f++, i, _T, *env);
      vectors.bits2_set(f++, i, _A, *env);
      vectors.bits2_set(f++, i, _A, *env);
      vectors.bits2_set(f++, i, _T, *env);
      vectors.bits2_set(f++, i, _T, *env);
      vectors.bits2_set(f++, i, _A, *env);
      vectors.bits2_set(f++, i, _T, *env);
      vectors.bits2_set(f++, i, UN, *env);
      vectors.bits2_set(f++, i, UN, *env);
      vectors.bits2_set(f++, i, _A, *env);
    }
    {
      const int _A = 0;
      const int _T = 1;
      int f = 0;
      const int i = 2;
      vectors.bits2_set(f++, i, _A, *env);
      vectors.bits2_set(f++, i, _A, *env);
      vectors.bits2_set(f++, i, _A, *env);
      vectors.bits2_set(f++, i, _T, *env);
      vectors.bits2_set(f++, i, _T, *env);
      vectors.bits2_set(f++, i, _A, *env);
      vectors.bits2_set(f++, i, UN, *env);
      vectors.bits2_set(f++, i, _A, *env);
      vectors.bits2_set(f++, i, _A, *env);
      vectors.bits2_set(f++, i, UN, *env);
    }
  }

  GMMetrics metrics_value = GMMetrics_null();
  GMMetrics* metrics = &metrics_value;
  MetricsMem metrics_mem(env_value);
  GMMetrics_create(metrics, env->data_type_metrics(), &dm,
                   &metrics_mem, env);

  if (env->is_proc_active())
    printf("%s\n", options.c_str());

  ComputeMetrics::compute(*metrics, vectors, *env);

  Checksum cksum;
  Checksum cksum_local;
  Checksum::compute(cksum, cksum_local, *metrics, *env);
  Driver::print_output_sync(cksum, *env);

  if (env->is_proc_active()) {
    const double result000 =
      Metrics_is_pass_threshold_noshrink(*metrics, 0, 0, 0, 0, *env) ?
      Metrics_ccc_duo_get_3(*metrics, 0, 0, 0, 0, *env) : 0;
    const double result001 =
      Metrics_is_pass_threshold_noshrink(*metrics, 0, 0, 0, 1, *env) ?
      Metrics_ccc_duo_get_3(*metrics, 0, 0, 0, 1, *env) : 0;
    const double result010 =
      Metrics_is_pass_threshold_noshrink(*metrics, 0, 0, 1, 0, *env) ?
      Metrics_ccc_duo_get_3(*metrics, 0, 0, 1, 0, *env) : 0;
    const double result011 =
      Metrics_is_pass_threshold_noshrink(*metrics, 0, 0, 1, 1, *env) ?
      Metrics_ccc_duo_get_3(*metrics, 0, 0, 1, 1, *env) : 0;

    const double result100 =
      Metrics_is_pass_threshold_noshrink(*metrics, 0, 1, 0, 0, *env) ?
      Metrics_ccc_duo_get_3(*metrics, 0, 1, 0, 0, *env) : 0;
    const double result101 =
      Metrics_is_pass_threshold_noshrink(*metrics, 0, 1, 0, 1, *env) ?
      Metrics_ccc_duo_get_3(*metrics, 0, 1, 0, 1, *env) : 0;
    const double result110 =
      Metrics_is_pass_threshold_noshrink(*metrics, 0, 1, 1, 0, *env) ?
      Metrics_ccc_duo_get_3(*metrics, 0, 1, 1, 0, *env) : 0;
    const double result111 =
      Metrics_is_pass_threshold_noshrink(*metrics, 0, 1, 1, 1, *env) ?
      Metrics_ccc_duo_get_3(*metrics, 0, 1, 1, 1, *env) : 0;

    //const double result000 = Metrics_ccc_duo_get_3(*metrics, 0, 0, 0, 0, *env);
    //const double result001 = Metrics_ccc_duo_get_3(*metrics, 0, 0, 0, 1, *env);
    //const double result010 = Metrics_ccc_duo_get_3(*metrics, 0, 0, 1, 0, *env);
    //const double result011 = Metrics_ccc_duo_get_3(*metrics, 0, 0, 1, 1, *env);
    //const double result100 = Metrics_ccc_duo_get_3(*metrics, 0, 1, 0, 0, *env);
    //const double result101 = Metrics_ccc_duo_get_3(*metrics, 0, 1, 0, 1, *env);
    //const double result110 = Metrics_ccc_duo_get_3(*metrics, 0, 1, 1, 0, *env);
    //const double result111 = Metrics_ccc_duo_get_3(*metrics, 0, 1, 1, 1, *env);

    printf("COMPUTED: A A A  %.17f\n", result000);
    printf("COMPUTED: A A T  %.17f\n", result001);
    printf("COMPUTED: A T A  %.17f\n", result010);
    printf("COMPUTED: A T T  %.17f\n", result011);
    printf("COMPUTED: T A A  %.17f\n", result100);
    printf("COMPUTED: T A T  %.17f\n", result101);
    printf("COMPUTED: T T A  %.17f\n", result110);
    printf("COMPUTED: T T T  %.17f\n", result111);
    printf("\n");

    const double s0_0 = 3;
    const double s0_1 = 4;
    const double c0 = s0_0 + s0_1;

    const double s1_0 = 4;
    const double s1_1 = 4;
    const double c1 = s1_0 + s1_1;

    const double s2_0 = 6;
    const double s2_1 = 2;
    const double c2 = s2_0 + s2_1;

    const double f0_0 = s0_0 / c0;
    const double f0_1 = s0_1 / c0;
    const double f1_0 = s1_0 / c1;
    const double f1_1 = s1_1 / c1;
    const double f2_0 = s2_0 / c2;
    const double f2_1 = s2_1 / c2;

    const double r0_000 = 0;
    const double r0_001 = 0;
    const double r0_010 = 1;
    const double r0_011 = 0;
    const double r0_100 = 0;
    const double r0_101 = 0;
    const double r0_110 = 0;
    const double r0_111 = 0;

    const double r1_000 = 0;
    const double r1_001 = 0;
    const double r1_010 = 0;
    const double r1_011 = 0;
    const double r1_100 = 1;
    const double r1_101 = 0;
    const double r1_110 = 0;
    const double r1_111 = 0;

    const double r2_000 = 0;
    const double r2_001 = 0;
    const double r2_010 = 0;
    const double r2_011 = 0;
    const double r2_100 = 1;
    const double r2_101 = 0;
    const double r2_110 = 0;
    const double r2_111 = 0;

    const double r4_000 = 0;
    const double r4_001 = 0;
    const double r4_010 = 0;
    const double r4_011 = 0;
    const double r4_100 = 0;
    const double r4_101 = 0;
    const double r4_110 = 0;
    const double r4_111 = 1;

    const double r5_000 = 0;
    const double r5_001 = 0;
    const double r5_010 = 0;
    const double r5_011 = 0;
    const double r5_100 = 1;
    const double r5_101 = 0;
    const double r5_110 = 0;
    const double r5_111 = 0;

    //const double r3_*** = 0
    //const double r6_*** = 0
    //const double r7_*** = 0
    //const double r8_*** = 0
    //const double r9_*** = 0

    const double r_000 = r0_000 + r1_000 + r2_000 + r4_000 + r5_000;
    const double r_001 = r0_001 + r1_001 + r2_001 + r4_001 + r5_001;
    const double r_010 = r0_010 + r1_010 + r2_010 + r4_010 + r5_010;
    const double r_011 = r0_011 + r1_011 + r2_011 + r4_011 + r5_011;
    const double r_100 = r0_100 + r1_100 + r2_100 + r4_100 + r5_100;
    const double r_101 = r0_101 + r1_101 + r2_101 + r4_101 + r5_101;
    const double r_110 = r0_110 + r1_110 + r2_110 + r4_110 + r5_110;
    const double r_111 = r0_111 + r1_111 + r2_111 + r4_111 + r5_111;

    const double c = r_000 + r_001 + r_010 + r_011 +
                     r_100 + r_101 + r_110 + r_111;

    const double f_000 = r_000 / c;
    const double f_001 = r_001 / c;
    const double f_010 = r_010 / c;
    const double f_011 = r_011 / c;
    const double f_100 = r_100 / c;
    const double f_101 = r_101 / c;
    const double f_110 = r_110 / c;
    const double f_111 = r_111 / c;

    const double fm = 4;
    const double cp = 2 / (double) 3;

    const double ref000 = fm * f_000 * (1-cp*f0_0) * (1-cp*f1_0) * (1-cp*f2_0);
    const double ref001 = fm * f_001 * (1-cp*f0_0) * (1-cp*f1_0) * (1-cp*f2_1);
    const double ref010 = fm * f_010 * (1-cp*f0_0) * (1-cp*f1_1) * (1-cp*f2_0);
    const double ref011 = fm * f_011 * (1-cp*f0_0) * (1-cp*f1_1) * (1-cp*f2_1);
    const double ref100 = fm * f_100 * (1-cp*f0_1) * (1-cp*f1_0) * (1-cp*f2_0);
    const double ref101 = fm * f_101 * (1-cp*f0_1) * (1-cp*f1_0) * (1-cp*f2_1);
    const double ref110 = fm * f_110 * (1-cp*f0_1) * (1-cp*f1_1) * (1-cp*f2_0);
    const double ref111 = fm * f_111 * (1-cp*f0_1) * (1-cp*f1_1) * (1-cp*f2_1);

    const double ref000t = ref000 > tLLL || (ref000 + ref111 > tLLLHHH &&
      ref000 > 0 && ref111 > 0) ? ref000 : (double)0;
    const double ref001t = ref001 > tLLH ? ref001 : (double)0;
    const double ref010t = ref010 > tLLH ? ref010 : (double)0;
    const double ref011t = ref011 > tLHH ? ref011 : (double)0;
    const double ref100t = ref100 > tLLH ? ref100 : (double)0;
    const double ref101t = ref101 > tLHH ? ref101 : (double)0;
    const double ref110t = ref110 > tLHH ? ref110 : (double)0;
    const double ref111t = ref111 > tHHH || (ref000 + ref111 > tLLLHHH &&
      ref000 > 0 && ref111 > 0) ? ref111 : (double)0;

    printf("EXPECTED: A A A  %.17f\n", ref000t);
    printf("EXPECTED: A A T  %.17f\n", ref001t);
    printf("EXPECTED: A T A  %.17f\n", ref010t);
    printf("EXPECTED: A T T  %.17f\n", ref011t);
    printf("EXPECTED: T A A  %.17f\n", ref100t);
    printf("EXPECTED: T A T  %.17f\n", ref101t);
    printf("EXPECTED: T T A  %.17f\n", ref110t);
    printf("EXPECTED: T T T  %.17f\n", ref111t);
    printf("\n");

    const double eps = 1.e-3;

    EXPECT_EQ(true, fabs(result000 - ref000t) < eps);
    EXPECT_EQ(true, fabs(result001 - ref001t) < eps);
    EXPECT_EQ(true, fabs(result010 - ref010t) < eps);
    EXPECT_EQ(true, fabs(result011 - ref011t) < eps);
    EXPECT_EQ(true, fabs(result100 - ref100t) < eps);
    EXPECT_EQ(true, fabs(result101 - ref101t) < eps);
    EXPECT_EQ(true, fabs(result110 - ref110t) < eps);
    EXPECT_EQ(true, fabs(result111 - ref111t) < eps);
  }

  GMMetrics_destroy(metrics, env);
  vectors.deallocate();
  dm.deallocate();
} // DriverTest_ccc3_simple_sparse_compute_method

//=============================================================================

void DriverTest_duo3_simple_sparse_(int compute_method) {

//  if (comet::ComputeMethod::REF == compute_method ||
//      comet::BuildHas::DOUBLE_PREC) {
//
//    // REF method and DP build do not remove threshold-failed values,
//    // so don't test that.
//    const double tLLL = .01;
//    const double tLLH = .01;
//    const double tLHH = .01;
//    const double tHHH = .01;
//    const double tLLLHHH = .01;
//
//    DriverTest_duo3_simple_sparse_(compute_method, tLLL, tLLH, tLHH, tHHH,
//      tLLLHHH);
//
//  } else {

    const double tLLL = .01;
    const double tLHH = .01;
    for (double tLLH : {.19, .20, .49, .50}) {
    for (double tHHH : {.27, .28}) {
    for (double tLLLHHH : {.27, .28}) {
      DriverTest_duo3_simple_sparse_(compute_method, tLLL, tLLH, tLHH, tHHH,
        tLLLHHH);
    }
    }
    }

// }
}

//=============================================================================

void DriverTest_duo3_simple_sparse_() {
  DriverTest_duo3_simple_sparse_(comet::ComputeMethod::REF);
  DriverTest_duo3_simple_sparse_(comet::ComputeMethod::CPU);
  DriverTest_duo3_simple_sparse_(comet::ComputeMethod::GPU);
}

//=============================================================================

void DriverTest_tc_() {

    char options1[1024];
    char options2[1024];

    char options_template[] =
        "--metric_type %s "
        "--num_proc_vector %i --num_field 100 --num_vector %i "
        //"--num_proc_vector %i --num_field 1 --num_vector %i "
        "--compute_method %s --sparse %s "
        "--problem_type random --verbosity %i --tc %i --num_way %i "
        "--num_tc_steps %i --all2all yes" ;

    const int num_proc_vector = comet::System::num_proc() >= 2 ? 2 : 1;
    //const int num_proc_vector = 1;

    typedef comet::TC TC;

    for (int is_duo=0; is_duo<=1; ++is_duo) {
    //for (int is_duo=1; is_duo<=1; ++is_duo) {
    for (int gpu=0; gpu<=1; ++gpu) {
    //for (int gpu=0; gpu<=0; ++gpu) {
    for (int num_tc_steps=1; num_tc_steps<=3; ++num_tc_steps) {
    //for (int num_tc_steps=1; num_tc_steps<=1; ++num_tc_steps) {
    for (int nv=10; nv<=10; ++nv) {
    //for (int nv=3; nv<=3; ++nv) {
    for (int num_way=2; num_way<=3; ++num_way) {
    //for (int num_way=3; num_way<=3; ++num_way) {
    for (int sparse=0; sparse<=1; ++sparse) {
    //for (int sparse=1; sparse<=1; ++sparse) {
      if (is_duo && 0 == sparse) continue;
    for (int tc=1; tc<TC::NUM; ++tc) {
    //for (int tc=4; tc<=4; ++tc) {
      if (nv/num_proc_vector < num_way) continue;
      if (!is_duo && TC::B1 == tc) continue;

      sprintf(options1, options_template, is_duo ? "duo" : "ccc",
              num_proc_vector, nv, "REF",
              sparse ? "yes" : "no", 1, 0, num_way, 1);
      sprintf(options2, options_template, is_duo ? "duo" : "ccc",
              num_proc_vector, nv, gpu ? "GPU" : "CPU",
              sparse ? "yes" : "no", 1, tc, num_way, num_tc_steps);
      if (is_duo && 3 == num_way && !is_using_tc(options2)) continue;
      EXPECT_EQ(true, compare_2runs(options1, options2));
    }
    }
    }
    }
    }
    }
    }
} // DriverTest_tc_

//=============================================================================

void DriverTest_subbyte_gemm_() {

  char options1[1024];
  char options2[1024];

  char options_template[] =
      "--metric_type duo --num_way %i --num_proc_vector 1 "
      "--num_vector %i --num_field %i "
      "--compute_method %s --sparse yes --all2all yes "
      "--problem_type random --verbosity %i --tc %i "
      "--num_tc_steps %i --threshold .05 --metrics_shrink 1.0 ";

  typedef comet::TC TC;

  for (int tc_method : {TC::INT4, TC::B1})
  for (int num_way : {2, 3})
  //for (int num_way : {2})
  for (int num_tc_steps : {1, 2})
  //for (int num_tc_steps : {1})
  for (int num_vector = num_way; num_vector < 5; ++num_vector) {
  //for (int num_vector = 2; num_vector < 3; ++num_vector) {
    //if (num_vector < num_way)
    //  continue;
  // Examine num_field values nearby possible block boundaries.
  int num_field_prev = 0;
  for (int nfbdry : {1, 2, 4, 8, 16, 32, 64, 128}) {
    const int range = 1;
  for (int nf = nfbdry - range; nf <= nfbdry + range; ++nf) {
  //for (int nf = nfbdry + 1; nf < nfbdry + 2; ++nf) {
    const int num_field = nf;
    // Skip if this num_field already visited.
    if (num_field <= num_field_prev)
      continue;

    //const int num_vector_ = 256; // 256;
    //const int num_field_ = 256;

    sprintf(options1, options_template, num_way,
            num_vector, num_field, "REF", 1, 0, 1);
    sprintf(options2, options_template, num_way,
            num_vector, num_field, "GPU", 1, tc_method, num_tc_steps);
    EXPECT_EQ(true, compare_2runs(options1, options2));
    num_field_prev = num_field;
  }
  }
  }

} // DriverTest_subbyte_gemm_

//=============================================================================

void DriverTest_threshold_() {

    char options1[1024];
    char options2[1024];

    const int num_proc_vector = comet::System::num_proc() >= 3 ? 3 : 1;
    //const int num_proc_vector = 1;

    // NOTE: --problem_type analytic can give high (nonphysical)
    // values of the metrics compared to --problem_type random.

    char options_template[] =
        "--metric_type %s "
        "--num_proc_vector %i "
        //"--num_field 71 --num_vector 17 " // takes too long for some cases.
        //"--num_field 23 --num_vector 13 "
        "--num_field 13 --num_vector 11 "
        //"--num_field 1 --num_vector 2 "
        //"--num_field 71 --num_vector 20 "
        "--num_way %i "
        "--all2all yes --sparse yes "
        //"--problem_type random "
        "--problem_type analytic "
        "--threshold %s "
        "--tc %i "
        "--compute_method %s "
        "--metrics_shrink %f "
        "--verbosity 1 ";

    typedef comet::ComputeMethod CM;
    typedef comet::MetricType MT;
    typedef comet::TC TC;

    // To help avoid problems with comparing floating point values.
    const double fuzz = .001;

    for (int thresholds_case : {0, 1})
    //for (int thresholds_case : {0})
    for (double metrics_shrink : {1., (10./3.)+fuzz})
    //XXXfor (int metrics_shrink_case : {0, 1}) {
    //XXX  const double metrics_shrink = 0 == metrics_shrink_case ? 1. :
    for (int num_way : {2, 3}) {
    //for (int num_way : {3}) {
      const double tmax = 2 == num_way ? .99 : .65;
    for (int metric_type : {MT::CCC, MT::DUO})
    //for (int metric_type : {MT::CCC})
    //for (double threshold : {tmax, tmax*.8, tmax*.5, .01, 0.})
    for (double threshold : {tmax, tmax*.8, .01, 0.})
    //for (double threshold : {tmax})
    for (int compute_method : {CM::CPU, CM::GPU})
    //for (int compute_method : {CM::GPU})
    for (int tc=1; tc<TC::NUM; ++tc) {
    //for (int tc=2; tc<=2; ++tc) {
      // Some cases have too many thresholded values to be able to shrink.
      const bool is_trying_shrink = metrics_shrink > 1.+fuzz;
      if (is_trying_shrink && 2 == num_way &&
          threshold < tmax-fuzz && threshold > 0+fuzz)
        continue;
      if (is_trying_shrink && (CM::GPU != compute_method ||
           threshold < .6 * tmax || comet::BuildHas::DOUBLE_PREC))
        continue;
      // Unimplemented.
      if (MT::CCC == metric_type && TC::B1 == tc) continue;

      // Some cases take very long to run.
      if (comet::BuildHas::DOUBLE_PREC && (1 != thresholds_case ||
         tmax*.8 != threshold || TC::AUTO == tc)) continue;
      //   (2 == tc || 4 == tc || 6 == tc))

      char thresholds1[1024];
      char thresholds2[1024];

      const double t = threshold;

      if (0 == thresholds_case) {
        sprintf(thresholds1, "%f", threshold);
        if (2 == num_way) {
          sprintf(thresholds2, "%f,%f,%f,1e6", t, t, t);
        } else {
          sprintf(thresholds2, "%f,%f,%f,%f,1e6", t, t, t, t);
        }
      } else {
        if (2 == num_way) {
          //sprintf(thresholds1, "%f,%f,%f,%f", t, t, t, t);
          //sprintf(thresholds2, "%f,%f,%f,%f", t, t, t, t);
          sprintf(thresholds1, "%f,%f,%f,%f", t, 1.05*t, 1.1*t, 2.1*t);
          sprintf(thresholds2, "%f,%f,%f,%f", t, 1.05*t, 1.1*t, 2.1*t);
        } else {
          //sprintf(thresholds1, "%f,%f,%f,%f,%f", t, t, t, t, t);
          //sprintf(thresholds2, "%f,%f,%f,%f,%f", t, t, t, t, t);
          sprintf(thresholds1, "%f,%f,%f,%f,%f", t, 1.05*t, 1.1*t, 1.15*t, 2.1*t);
          sprintf(thresholds2, "%f,%f,%f,%f,%f", t, 1.05*t, 1.1*t, 1.15*t, 2.1*t);
        }
      }

      sprintf(options1, options_template, MT::str(metric_type),
        num_proc_vector, num_way, thresholds1, TC::NO, CM::str(CM::REF), 1.);
      sprintf(options2, options_template, MT::str(metric_type),
        num_proc_vector, num_way, thresholds2, tc, CM::str(compute_method),
        metrics_shrink);
      test_2runs(options1, options2);
    } 
    }
} // DriverTest_tc_

//=============================================================================

void DriverTest_file_output_() {

    char options1[1024];
    char options2[1024];

    //const int num_proc_vector = 5;
    //const int num_proc_vector = 2;
    const int num_proc_vector = comet::System::num_proc() >= 2 ? 11 : 1;

    char options_template[] =
        "--metric_type %s "
        "--num_proc_vector %i "
        "--num_field 19 --num_vector 55 "
        "--num_way %i "
        "--all2all %s "
        "--sparse %s "
        "--problem_type random "
        "--threshold .60 "
        "--tc %i "
        "--compute_method %s "
        "--num_phase %i "
        "--num_stage %i "
        "--verbosity 1 "
        "%s";

    typedef comet::ComputeMethod CM;
    typedef comet::MetricType MT;
    typedef comet::TC TC;
    typedef comet::NumWay NumWay;

    // Metrics output file.

    for (int num_way : {NumWay::_2, NumWay::_3})
    for (int all2all : {1})
    for (int metric_type : {MT::CZEK, MT::CCC, MT::DUO})
    for (int compute_method : {CM::CPU, CM::GPU}) {


      // Some cases take very long to run.
      if (comet::BuildHas::DOUBLE_PREC && NumWay::_3 == num_way &&
          COMET_COMPUTE_CAPABILITY >= 800) continue;

      const int num_stage = NumWay::_3 == num_way && all2all ? 2 : 1;
      const int num_phase = 1 == num_proc_vector ? 1 : all2all ? 2 : 1;
      const char* sparse = MT::CZEK == metric_type ? "no" : "yes";

      sprintf(options1, options_template,
        MT::str(metric_type), num_proc_vector, num_way, all2all ? "yes" : "no",
        sparse, TC::NO, CM::str(CM::REF), num_phase, num_stage, "");

      sprintf(options2, options_template,
        MT::str(metric_type), num_proc_vector, num_way, all2all ? "yes" : "no",
        sparse, TC::AUTO, CM::str(compute_method), num_phase, num_stage,
        "--output_file_stub DriverTest_file_output");

      test_2runs(options1, options2);
    }

    // Histograms.

    for (int num_way : {NumWay::_2, NumWay::_3})
    for (int all2all : {1})
    for (int metric_type : {MT::CCC, MT::DUO})
    //for (int metric_type : {MT::DUO})
    for (int compute_method : {CM::GPU})
    for (int tc : {TC::NO, TC::AUTO}) {
      // Some cases take very long to run.
      if (comet::BuildHas::DOUBLE_PREC && NumWay::_3 == num_way &&
          COMET_COMPUTE_CAPABILITY >= 800) continue;

      const int num_stage = NumWay::_3 == num_way && all2all ? 2 : 1;
      const int num_phase = 1 == num_proc_vector ? 1 : all2all ? 2 : 1;
      const char* sparse = MT::CZEK == metric_type ? "no" : "yes";

      char hfile_template[] =
          "--histograms_file DriverTest_file_output_histogram_%i.csv";

      char hfile[1024];

      sprintf(hfile, hfile_template, num_way);

      sprintf(options1, options_template,
        MT::str(metric_type), num_proc_vector, num_way, all2all ? "yes" : "no",
        sparse, TC::NO, CM::str(CM::REF), num_phase, num_stage, "");

      sprintf(options2, options_template,
        MT::str(metric_type), num_proc_vector, num_way, all2all ? "yes" : "no",
        sparse, tc, CM::str(compute_method), num_phase, num_stage,
        hfile);
        //"--histograms_file DriverTest_file_output_histogram.csv");

      test_2runs(options1, options2);
    }

} // DriverTest_file_output_

//=============================================================================

void DriverTest_ccc2_duo2_(const char* const metric_type) {
  char options1[1024];
  char options2[1024];
  char options3[1024];
#if 1

  //----------
  //---2-way, all2all no
  //----------

  char options_template_1[] =
      "--sparse yes "
      "--metric_type %s --verbosity %i "
      //"--problem_type analytic " 
      "--num_proc_vector 1 --num_field %i --num_vector_local %i "
      "--compute_method %s";

  sprintf(options1, options_template_1, metric_type, 2, 1, 2, "REF");
  sprintf(options2, options_template_1, metric_type, 2, 1, 2, "CPU");
  sprintf(options3, options_template_1, metric_type, 2, 1, 2, "GPU");
  EXPECT_EQ(true, compare_3runs(options1, options2, options3));

  sprintf(options1, options_template_1, metric_type, 1, 100, 2, "REF");
  sprintf(options2, options_template_1, metric_type, 1, 100, 2, "CPU");
  sprintf(options3, options_template_1, metric_type, 1, 100, 2, "GPU");
  EXPECT_EQ(true, compare_3runs(options1, options2, options3));

  sprintf(options1, options_template_1, metric_type, 1, 100, 48, "REF");
  sprintf(options2, options_template_1, metric_type, 1, 100, 48, "CPU");
  sprintf(options3, options_template_1, metric_type, 1, 100, 48, "GPU");
  EXPECT_EQ(true, compare_3runs(options1, options2, options3));

  for (int i = 1; i <= 100; ++i) {
    if ((i+3) % 32 > 6) {
      continue;
    }
    sprintf(options1, options_template_1, metric_type, 1, i, 48, "REF");
    sprintf(options2, options_template_1, metric_type, 1, i, 48, "CPU");
    sprintf(options3, options_template_1, metric_type, 1, i, 48, "GPU");
    EXPECT_EQ(true, compare_3runs(options1, options2, options3));
  }

  //----------
  //---2-way, all2all yes, small
  //----------

  char options_template_2[] =
      "--sparse yes "
      "--metric_type %s --verbosity %i "
      //"--problem_type analytic "
      "--num_proc_vector %i --num_field %i --num_vector_local %i "
      "--compute_method %s --all2all %s";

  sprintf(options1, options_template_2, metric_type, 2, 1, 1, 2,
          "REF", "no");
  sprintf(options2, options_template_2, metric_type, 2, 1, 1, 2,
          "CPU", "yes");
  sprintf(options3, options_template_2, metric_type, 2, 1, 1, 2,
          "GPU", "yes");
  EXPECT_EQ(true, compare_3runs(options1, options2, options3));

  sprintf(options1, options_template_2, metric_type, 1, 1, 1, 4,
          "REF", "no");
  sprintf(options2, options_template_2, metric_type, 1, 2, 1, 2,
          "CPU", "yes");
  sprintf(options3, options_template_2, metric_type, 1, 2, 1, 2,
          "GPU", "yes");
  EXPECT_EQ(true, compare_3runs(options1, options2, options3));

  char options_template_2a[] =
      "--sparse yes "
      "--num_proc_vector 1 --num_field 2 --num_vector 5 "
      //"--problem_type analytic "
      "--compute_method REF --all2all yes --metric_type %s";

  char options_template_2b[] =
      "--sparse yes "
      "--num_proc_vector 2 --num_field 2 --num_vector 5 "
      //"--problem_type analytic "
      "--compute_method CPU --all2all yes --metric_type %s";

  char options_template_2c[] =
      "--sparse yes "
      "--num_proc_vector 2 --num_field 2 --num_vector 5 "
      //"--problem_type analytic "
      "--compute_method GPU --all2all yes --metric_type %s";

  sprintf(options1, options_template_2a, metric_type);
  sprintf(options2, options_template_2b, metric_type);
  sprintf(options3, options_template_2c, metric_type);
  EXPECT_EQ(true, compare_3runs(options1, options2, options3));

  char options_template_2a2[] =
      "--sparse yes "
      "--num_proc_vector 1 --num_proc_field 1 "
      "--num_field 7 --num_vector 2 "
      //"--problem_type analytic "
      "--compute_method REF --all2all yes --metric_type %s";

  char options_template_2b2[] =
      "--sparse yes "
      "--num_proc_vector 1 --num_proc_field 3 "
      "--num_field 7 --num_vector 2 "
      //"--problem_type analytic "
      "--compute_method GPU --all2all yes --metric_type %s";

  sprintf(options1, options_template_2a2, metric_type);
  sprintf(options2, options_template_2b2, metric_type);
  EXPECT_EQ(true, compare_2runs(options1, options2));

  //----------
  //---2-way, all2all yes, large
  //----------

  sprintf(options1, options_template_2, metric_type, 1, 1, 100, 48,
          "REF", "no");
  sprintf(options2, options_template_2, metric_type, 1, 1, 100, 48,
          "CPU", "yes");
  sprintf(options3, options_template_2, metric_type, 1, 1, 100, 48,
          "GPU", "yes");
  EXPECT_EQ(true, compare_3runs(options1, options2, options3));

  sprintf(options1, options_template_2, metric_type, 1, 1, 100, 48,
          "REF", "no");
  sprintf(options2, options_template_2, metric_type, 1, 2, 100, 24,
          "CPU", "yes");
  sprintf(options3, options_template_2, metric_type, 1, 2, 100, 24,
          "GPU", "yes");
  EXPECT_EQ(true, compare_3runs(options1, options2, options3));

  //----------
  //---num_proc_field
  //----------

  char options_template_5[] =
      "--sparse yes "
      "--metric_type %s --verbosity %i "
      "--num_proc_vector %i --num_proc_field %i "
      " --num_field %i --num_vector_local %i "
      //"--problem_type analytic "
      "--compute_method %s --all2all %s";

  for (int i = 1; i <= 3; ++i) {
    sprintf(options1, options_template_5, metric_type, 1, 1, 1, 60, 48, "REF", "no");
    sprintf(options2, options_template_5, metric_type, 1, 1, 2 * i - 1, 60, 48, "GPU",
            "yes");
    sprintf(options3, options_template_5, metric_type, 1, 1, 2 * i, 60, 48, "GPU", "yes");
    EXPECT_EQ(true, compare_3runs(options1, options2, options3));
    sprintf(options1, options_template_5, metric_type, 1, 1, 1, 60, 48, "REF", "no");
    sprintf(options2, options_template_5, metric_type, 1, 1, i, 60, 48, "GPU", "yes");
    sprintf(options3, options_template_5, metric_type, 1, 2, i, 60, 24, "GPU", "yes");
    EXPECT_EQ(true, compare_3runs(options1, options2, options3));
  }

  //----------
  //---num_proc_repl, 2-way
  //----------

  char options_template_10[] =
      "--sparse yes "
      "--metric_type %s "
      //"--problem_type analytic "
      "--num_field 4 --num_vector_local %i --compute_method %s --all2all yes "
      "--num_proc_vector %i --num_proc_repl %i "
      "--num_proc_field %i --num_way %i --num_stage %i";

  int num_vector_local = 0;
  int num_proc_vector = 0;
  int num_proc_repl = 0;
  int gpu = 0;
  int num_stage = 1;

  for (gpu=0; gpu<=1; ++gpu) {
    for (num_vector_local=4; num_vector_local<=5; ++num_vector_local) {
      for (num_proc_vector=1; num_proc_vector<=4; ++num_proc_vector) {
        for (num_proc_repl=1; num_proc_repl<=5; ++num_proc_repl) {
          const int num_proc_field = gpu ? 2 : 1;
          const int num_way = 2;
          sprintf(options1, options_template_10, metric_type,
                  num_vector_local*num_proc_vector,
                  "GPU", 1, 1,
                  1, num_way, 1);
          sprintf(options2, options_template_10, metric_type, num_vector_local,
                  gpu ? "GPU" : "CPU", num_proc_vector, num_proc_repl,
                  num_proc_field, num_way, num_stage);
          test_2runs(options1, options2);
        }
      }
    }
  }

  //----------
  //---ccc_param
  //----------

  char options_template_11a[] =
      "--sparse yes "
      "--metric_type %s --verbosity %i "
      //"--problem_type analytic "
      "--num_proc_vector 1 --num_field 30 --num_vector_local 40 "
      "--compute_method GPU";
  char options_template_11b[] =
      "--sparse yes "
      "--metric_type %s --verbosity %i "
      //"--problem_type analytic "
      "--num_proc_vector 1 --num_field 30 --num_vector_local 40 "
      "--compute_method GPU --ccc_param %.20e --problem_type analytic";

  sprintf(options1, options_template_11a, metric_type, 1);
  sprintf(options2, options_template_11b, metric_type, 1, ((double)2) / ((double)3));
  EXPECT_EQ(true, compare_2runs(options1, options2));

  int proc_num = 0;
  COMET_MPI_SAFE_CALL(MPI_Comm_rank(MPI_COMM_WORLD, &proc_num));

  sprintf(options1, options_template_11a, metric_type, 1);
  sprintf(options2, options_template_11b, metric_type, 1, ((double)1) / ((double)2));
  if (can_run(options1) && can_run(options2)) {
    const int result11 = compare_2runs(options1, options2);
    EXPECT_EQ(true, 0==proc_num ? ! result11 : true);
  }
#endif
#if 1

  //----------
  //---num_phase, 2-way
  //----------

  int num_phase = 0;

  for (int num_proc_vector=1; num_proc_vector<=8; ++num_proc_vector) {
    for (int num_proc_repl=1; num_proc_repl<=8; ++num_proc_repl) {
      for (int num_phase=2; num_phase<=8; ++num_phase) {
        if (!(num_phase <= 1 + num_proc_vector/2)) {
          continue;
        }
        char options_template[] =
          "--sparse yes "
          "--metric_type %s "
          //"--problem_type analytic "
          "--num_field 7 --num_vector 12 --compute_method GPU --all2all yes "
          "--num_proc_vector %i --num_proc_repl %i --num_phase %i "
          "--num_way 2";
        sprintf(options1, options_template, metric_type, num_proc_vector, num_proc_repl,
                num_phase);
        sprintf(options2, options_template, metric_type, 1, 1, 1);
        test_2runs(options1, options2);
      }
    }
  }

  //----------
  //---2-way, all2all yes, large, sparse
  //----------

  {
    char options_template_2[] =
      "--sparse yes "
      "--metric_type %s --verbosity %i "
      //"--problem_type analytic "
      "--num_proc_vector %i --num_field %i --num_vector %i "
      "--compute_method %s --all2all %s";

    const int nf = 100;
    const int nv = 48;

    sprintf(options1, options_template_2, metric_type, 1, 1, nf, nv, "REF", "no");
    sprintf(options2, options_template_2, metric_type, 1, 1, nf, nv, "CPU", "yes");
    sprintf(options3, options_template_2, metric_type, 1, 1, nf, nv, "GPU", "yes");
    EXPECT_EQ(true, compare_3runs(options1, options2, options3));

    sprintf(options1, options_template_2, metric_type, 1, 1, nf, nv, "REF", "no");
    sprintf(options2, options_template_2, metric_type, 1, 2, nf, nv, "CPU", "yes");
    sprintf(options3, options_template_2, metric_type, 1, 2, nf, nv, "GPU", "yes");
    EXPECT_EQ(true, compare_3runs(options1, options2, options3));
  }

  //----------
  //---file input
  //----------

  for (int num_vector = 13; num_vector <= 13; ++num_vector) {
    for (int num_field = 1; num_field <= 3*300; ++num_field) {
      create_vectors_file("ccc_duo_2way_in.bin", num_field, num_vector,
                          comet::MetricType::CCC, 2, comet::ProblemType::DEFAULT, 1);
      for (int num_proc_vector=3; num_proc_vector<=3; ++num_proc_vector) {
        for (int num_proc_field=1; num_proc_field<=3; ++num_proc_field) {

          const bool test_only_critical_values = true;

          if (test_only_critical_values) {
            const bool is_nearly_multiple_of_32 =
              (num_field/num_proc_field+4) % 32 <= 8;
            if (! is_nearly_multiple_of_32) {
              continue;
            }
          }

          char options1[1024];
          char options2[1024];

          char options_template[] =
                 "--sparse yes "
                 "--metric_type %s "
                 "--num_vector %i --num_field %i "
                 "--num_proc_vector %i --num_proc_field %i "
                 "--compute_method GPU --all2all yes %s --verbosity 1";

          sprintf(options1, options_template, metric_type,
                  num_vector, num_field, num_proc_vector, num_proc_field, "");

          sprintf(options2, options_template, metric_type,
                  num_vector, num_field, num_proc_vector, num_proc_field,
                  "--input_file ccc_duo_2way_in.bin");

          test_2runs(options1, options2);
        }
      }
    }
  }
#endif

} // DriverTest_ccc2_duo2_

//=============================================================================

void DriverTest_ccc3_duo3_(const char* const metric_type) {
  char options1[1024];
  char options2[1024];
#if 1
  char options3[1024];
//  const bool is_duo = 'd' == metric_type[0];

  //----------
  //---3-way, all2all no
  //----------

  char options_template_3[] =
      "--sparse yes "
      "--metric_type %s --verbosity %i "
      "--num_proc_vector 1 --num_field %i --num_vector_local %i "
      //"--problem_type analytic "
      "--compute_method %s --num_way 3";

//  if (!is_duo) {
    sprintf(options1, options_template_3, metric_type, 2, 1, 3, "REF");
    sprintf(options2, options_template_3, metric_type, 2, 1, 3, "CPU");
    sprintf(options3, options_template_3, metric_type, 2, 1, 3, "GPU");
    EXPECT_EQ(true, compare_3runs(options1, options2, options3));

    sprintf(options1, options_template_3, metric_type, 1, 100, 48, "REF");
    sprintf(options2, options_template_3, metric_type, 1, 100, 48, "CPU");
    sprintf(options3, options_template_3, metric_type, 1, 100, 48, "GPU");
    EXPECT_EQ(true, compare_3runs(options1, options2, options3));

    for (int i = 1; i <= 100; ++i) {
      if ((i+3) % 32 > 6) {
        continue;
      }
      sprintf(options1, options_template_3, metric_type, 0, i, 24, "REF");
      sprintf(options2, options_template_3, metric_type, 0, i, 24, "CPU");
      sprintf(options3, options_template_3, metric_type, 0, i, 24, "GPU");
      EXPECT_EQ(true, compare_3runs(options1, options2, options3));
    }
//  }

  //----------
  //---3-way, all2all yes, small
  //----------

  char options_template_4[] =
      "--sparse yes "
      "--metric_type %s --verbosity %i "
      "--num_proc_vector %i --num_field %i --num_vector_local %i "
      //"--problem_type analytic "
      "--compute_method %s --all2all %s --num_way 3";

//  if (!is_duo) {
    sprintf(options1, options_template_4, metric_type, 2, 1, 1, 3, "REF", "no");
    sprintf(options2, options_template_4, metric_type, 2, 1, 1, 3, "CPU", "yes");
    sprintf(options3, options_template_4, metric_type, 2, 1, 1, 3, "GPU", "yes");
    EXPECT_EQ(true, compare_3runs(options1, options2, options3));

    sprintf(options1, options_template_4, metric_type, 1, 1, 1, 6, "REF", "no");
    sprintf(options2, options_template_4, metric_type, 1, 2, 1, 3, "CPU", "yes");
    sprintf(options3, options_template_4, metric_type, 1, 2, 1, 3, "GPU", "yes");
    EXPECT_EQ(true, compare_3runs(options1, options2, options3));
//  }

  //----------
  //---3-way, all2all yes, large
  //----------

//  if (!is_duo) {
    sprintf(options1, options_template_4, metric_type, 1, 1, 100, 48, "REF", "no");
    sprintf(options2, options_template_4, metric_type, 1, 1, 100, 48, "CPU", "yes");
    sprintf(options3, options_template_4, metric_type, 1, 1, 100, 48, "GPU", "yes");
    EXPECT_EQ(true, compare_3runs(options1, options2, options3));

    sprintf(options1, options_template_4, metric_type, 1, 1, 100, 48, "REF", "no");
    sprintf(options2, options_template_4, metric_type, 1, 4, 100, 12, "CPU", "yes");
    sprintf(options3, options_template_4, metric_type, 1, 4, 100, 12, "GPU", "yes");
    EXPECT_EQ(true, compare_3runs(options1, options2, options3));

    char options_template_4a[] =
        "--sparse yes "
        "--num_proc_vector 1 --num_field 2 --num_vector 16 "
        //"--problem_type analytic "
        "--compute_method REF --all2all yes --num_way 3 "
        "--metric_type %s";

    char options_template_4b[] =
        "--sparse yes "
        "--num_proc_vector 3 --num_field 2 --num_vector 16 "
        //"--problem_type analytic "
        "--compute_method CPU --all2all yes --num_way 3 "
        "--metric_type %s";

    char options_template_4c[] =
        "--sparse yes "
        "--num_proc_vector 3 --num_field 2 --num_vector 16 "
        //"--problem_type analytic "
        "--compute_method GPU --all2all yes --num_way 3 "
        "--metric_type %s";

    sprintf(options1, options_template_4a, metric_type);
    sprintf(options2, options_template_4b, metric_type);
    sprintf(options3, options_template_4c, metric_type);
    EXPECT_EQ(true, compare_3runs(options1, options2, options3));

    char options_template_4a2[] =
        "--sparse yes "
        "--num_proc_vector 1 --num_proc_field 1 "
        //"--problem_type analytic "
        "--num_field 13 --num_vector 3 "
        "--compute_method REF --all2all yes --num_way 3 "
        "--metric_type %s";

    char options_template_4b2[] =
        "--sparse yes "
        "--num_proc_vector 1 --num_proc_field 5 "
        //"--problem_type analytic "
        "--num_field 13 --num_vector 3 "
        "--compute_method GPU --all2all yes --num_way 3 "
        "--metric_type %s";

    sprintf(options1, options_template_4a2, metric_type);
    sprintf(options2, options_template_4b2, metric_type);
    EXPECT_EQ(true, compare_2runs(options1, options2));
//  }

  //----------
  //---num_proc_repl, num_stage, 3-way
  //----------

  char options_template_10[] =
      "--sparse yes "
      "--metric_type %s "
      //"--problem_type analytic "
      "--num_field 4 --num_vector_local %i --compute_method %s --all2all yes "
      "--num_proc_vector %i --num_proc_repl %i "
      "--num_proc_field %i --num_way %i --num_stage %i";

//  if (!is_duo) {
    for (int gpu=0; gpu<=1; ++gpu) {
      for (int num_vector_local=6; num_vector_local<=18; num_vector_local+=12) {
        for (int num_proc_vector=1; num_proc_vector<=4; ++num_proc_vector) {
          for (int num_proc_repl=1; num_proc_repl<=5; ++num_proc_repl) {
            for (int num_stage=1; num_stage<=6; num_stage+=4) {
              const int num_proc_field = gpu ? 2 : 1;
              const int num_way = 3;
              sprintf(options1, options_template_10, metric_type,
                      num_vector_local*num_proc_vector,
                      "GPU", 1, 1,
                      1, num_way, 1);
              sprintf(options2, options_template_10, metric_type,
                      num_vector_local,
                      gpu ? "GPU" : "CPU", num_proc_vector, num_proc_repl,
                      num_proc_field, num_way, num_stage);
              test_2runs(options1, options2);
            }
          }
        }
      }
    }
//  }

  //----------
  //---3-way, all2all yes, large, sparse
  //----------

//  if (!is_duo) {
    char options_template_5[] =
        "--sparse yes "
        "--metric_type %s --verbosity %i "
        //"--problem_type analytic "
        "--num_proc_vector %i --num_field %i --num_vector_local %i "
        "--compute_method %s --all2all %s --num_way 3";

    sprintf(options1, options_template_5, metric_type, 1, 1, 100, 48, "REF", "no");
    sprintf(options2, options_template_5, metric_type, 1, 1, 100, 48, "CPU", "yes");
    sprintf(options3, options_template_5, metric_type, 1, 1, 100, 48, "GPU", "yes");
    EXPECT_EQ(true, compare_3runs(options1, options2, options3));

    sprintf(options1, options_template_5, metric_type, 1, 1, 100, 48, "REF", "no");
    sprintf(options2, options_template_5, metric_type, 1, 4, 100, 12, "CPU", "yes");
    sprintf(options3, options_template_5, metric_type, 1, 4, 100, 12, "GPU", "yes");
    EXPECT_EQ(true, compare_3runs(options1, options2, options3));
//  }
#endif

} // DriverTest_ccc3_duo3_

//=============================================================================

void DriverTest_ccc2_() {
  DriverTest_ccc2_duo2_("ccc");
}

//=============================================================================

void DriverTest_duo2_() {
  DriverTest_ccc2_duo2_("duo");
}

//=============================================================================

void DriverTest_ccc3_() {
  DriverTest_ccc3_duo3_("ccc");
}

//=============================================================================

void DriverTest_duo3_() {
  DriverTest_ccc3_duo3_("duo");
}

//=============================================================================

BEGIN_TESTS

#if 1
TEST(DriverTest, czek2) {
  DriverTest_czek2_();
}

TEST(DriverTest, czek3) {
  DriverTest_czek3_();
}

TEST(DriverTest, threshold) {
  DriverTest_threshold_();
}
#endif

#if 1
TEST(DriverTest, file_output) {
  DriverTest_file_output_();
}

TEST(DriverTest, subbyte_gemm) {
  DriverTest_subbyte_gemm_();
}

TEST(DriverTest, tc) {
  DriverTest_tc_();
}

TEST(DriverTest, ccc3_simple) {
  DriverTest_ccc3_simple_();
}

TEST(DriverTest, ccc3_simple_sparse) {
  DriverTest_ccc3_simple_sparse_();
}

TEST(DriverTest, duo3_simple_sparse) {
  DriverTest_duo3_simple_sparse_();
}

TEST(DriverTest, ccc2_simple) {
  DriverTest_ccc2_simple_();
}

TEST(DriverTest, ccc2_simple_sparse) {
  DriverTest_ccc2_simple_sparse_();
}

TEST(DriverTest, duo2_simple_sparse) {
  DriverTest_duo2_simple_sparse_();
}
#endif

#if 1
TEST(DriverTest, ccc2) {
  DriverTest_ccc2_();
}

TEST(DriverTest, ccc3) {
  DriverTest_ccc3_();
}

TEST(DriverTest, duo2) {
  DriverTest_duo2_();
}

TEST(DriverTest, duo3) {
  DriverTest_duo3_();
}
#endif

END_TESTS

//=============================================================================

GTEST_API_ int main(int argc, char** argv) {

# ifdef COMET_USE_GTEST
    ::testing::InitGoogleTest(&argc, argv);
# endif

  COMET_MPI_SAFE_CALL(MPI_Init(&argc, &argv));

  int comm_rank = 0;
  COMET_MPI_SAFE_CALL(MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank));

  if (comm_rank != 0) {
#   ifdef COMET_USE_GTEST
      ::testing::TestEventListeners& listeners =
        ::testing::UnitTest::GetInstance()->listeners();
      delete listeners.Release(listeners.default_result_printer());
#   endif
  }

  int result = RUN_ALL_TESTS();

  int result_g = 11;

  COMET_MPI_SAFE_CALL(MPI_Allreduce(&result, &result_g, 1, MPI_INT, MPI_MAX,
    MPI_COMM_WORLD));

  COMET_MPI_SAFE_CALL(MPI_Finalize());
  return result_g;
}

//=============================================================================

//-----------------------------------------------------------------------------
