//-----------------------------------------------------------------------------
/*!
 * \file   test_problems.hh
 * \author Wayne Joubert
 * \date   Mon Aug  7 17:02:51 EDT 2017
 * \brief  Generator for synthetic test problems, header.
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

#ifndef _COMET_TEST_PROBLEMS_HH_
#define _COMET_TEST_PROBLEMS_HH_

#include "env.hh"
#include "vectors.hh"
#include "metrics.hh"

#include "driver.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------

struct ProblemType {
  enum {
    RANDOM = 1,
    ANALYTIC = 2,
    DEFAULT = ANALYTIC
  };
};

//-----------------------------------------------------------------------------

class TestProblem {
public:
  static void set_vectors_synthetic(GMVectors* vectors, int problem_type,
                                    int verbosity, CEnv* env);

  static void check_metrics(GMMetrics* metrics, Driver& driver, CEnv* env);
};

//-----------------------------------------------------------------------------
/// \brief Help[er class with test problem size info.

struct TestProblemInfo {

  enum {NUM_SHUFFLE = 3};

  TestProblemInfo(size_t nva, size_t nfa, CEnv& env)
    : nva_(nva)
    , nfa_(nfa)
      // Upper bound on integer representable exactly by floating point type.
      // Account for cast to float in magma Volta version.
    , max_float_(((size_t)1) <<  
                 (env.data_type_vectors() == GM_DATA_TYPE_FLOAT ?
                 mantissa_digits<float>() : mantissa_digits<GMFloat>()))
      // Czek account for number of terms summed in denom or num
    , overflow_limit_(env.data_type_vectors() != GM_DATA_TYPE_FLOAT ? 1 :
                      env.num_way() == NumWay::_2 ? 2 : 4)
      // Sum nfa times down the vector, is it still exact.
    , value_limit_((max_float_ - 1) / (overflow_limit_ * nfa_))
    , value_min_(1)
    , value_max_(utils::min(value_min_+nva_, value_limit_))
    , num_group_(1 << NUM_SHUFFLE)
    , group_size_max_(utils::ceil(nfa_, num_group_))
  {
  }

  const size_t nva_;
  const size_t nfa_;
  const size_t max_float_;
  const size_t overflow_limit_;
  const size_t value_limit_;
  const size_t value_min_;
  const size_t value_max_;
  const size_t num_group_;
  const size_t group_size_max_;
};

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_TEST_PROBLEMS_HH_

//-----------------------------------------------------------------------------
