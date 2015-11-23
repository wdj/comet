/*---------------------------------------------------------------------------*/
/*!
 * \file   system_test.cc
 * \author Wayne Joubert
 * \date   Fri Nov  6 18:18:21 EST 2015
 * \brief  Perform high-level system tests.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include <stdio.h>

#include "gtest/gtest.h"

#include "env.h"
#include "vectors.h"
#include "metrics.h"
#include "compute_metrics.h"

/*===========================================================================*/
/*---Set the entries of the vectors---*/

void input_vectors(GMVectors* vectors, GMEnv* env) {
  switch (Env_data_type(env)) {
    /*--------------------*/
    case GM_DATA_TYPE_FLOAT: {
      int vector_local;
      for (vector_local = 0; vector_local < vectors->num_vector_local;
           ++vector_local) {
        int field;
        for (field = 0; field < vectors->num_field; ++field) {
          /*---compute element unique id---*/
          const size_t uid = field +
                         vectors->num_field * (vector_local +
                                               vectors->num_vector_local *
                                                  ((size_t)Env_proc_num(env)));
          /*---Generate large random number---*/
          size_t rand1 = uid;
          rand1 = gm_randomize(rand1);
          rand1 = gm_randomize(rand1);
          size_t rand2 = uid;
          rand2 = gm_randomize(rand2);
          rand2 = gm_randomize(rand2);
          rand2 = gm_randomize(rand2);
          size_t rand_value = rand1 + gm_randomize_max() * rand2;
          /*---Reduce so that after summing num_field times the integer
               still fully fits in double precision fraction part---*/
          rand_value >>= (64-52) + gm_log2(vectors->num_field);
          /*---Store as floating point value---*/
          GMFloat value = rand_value;
          GMVectors_float_set(vectors, field, vector_local, value, env);
        } /*---field---*/
      }   /*---vector_local---*/
    } break;
    /*--------------------*/
    case GM_DATA_TYPE_BIT: {
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

      int vector_local;
      for (vector_local = 0; vector_local < vectors->num_vector_local;
           ++vector_local) {
        int field;
        for (field = 0; field < vectors->num_field; ++field) {
          /*---compute element unique id---*/
          size_t index = field +
                         vectors->num_field * (vector_local +
                                               vectors->num_vector_local *
                                                  ((size_t)Env_proc_num(env)));
          /*---randomize---*/
          index = gm_randomize(index);
          /*---Calculate random number between 0 and 1---*/
          GMFloat rand_value = index / (GMFloat)gm_randomize_max();
          /*---Create single bit value---*/
          _Bool value = rand_value < .5 ? GM_BOOL_FALSE : GM_BOOL_TRUE;
          GMVectors_bit_set(vectors, field, vector_local, value, env);
        } /*---field---*/
      }   /*---vector_local---*/

    } break;
    /*--------------------*/
    default:
      GMAssert(GM_BOOL_FALSE ? "Invalid data type." : 0);
  }
}

/*===========================================================================*/

TEST(SystemTest,One) {

  GMEnv env = GMEnv_null();
  GMEnv_create(&env);

  if (Env_proc_num(&env) == 0) {
    printf("///////// Running main()  /////////\n");
  }

  GMEnv_destroy(&env);

  EXPECT_EQ(1, 1);
}

/*===========================================================================*/

GTEST_API_ int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  testing::InitGoogleTest(&argc, argv);

  int result = RUN_ALL_TESTS();

  MPI_Finalize();
  return result;
}

/*===========================================================================*/

/*---------------------------------------------------------------------------*/
