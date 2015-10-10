/****************************************************************/
/* 3-Way Czekanowski Similarity Metric                          */
/*                                                              */
/* Code Author: James Nance                                     */
/* Created: October, 2015                                       */
/****************************************************************/

#include <stdio.h> /*FIX*/
#include <stdlib.h>

#include "magma_minproduct.h"
#include "magma_minproduct_lapack.h"

#include "env.h"
#include "vectors.h"
#include "metrics.h"
#include "compute_metrics_czekanowski_3way.h"

/*===========================================================================*/

static void compute_vector_sums(Vectors* vectors,
                                Float_t* __restrict__ vector_sums,
                                Env* env) {
  Assert(vectors != NULL);
  Assert(vector_sums != NULL);
  Assert(env != NULL);

  int i = 0;
  for ( i = 0; i < vectors->num_vector_local; ++i ) {
    Float_t sum = 0;
    int field = 0;
    for ( field = 0; field < vectors->num_field; ++field ) {
      Float_t value = Vectors_float_get(vectors, field, i, env);
      sum += value;
    }
    vector_sums[i] = sum;
  }
}

/*===========================================================================*/

void compute_metrics_czekanowski_3way_cpu(Metrics* metrics,
                                          Vectors* vectors,
                                          Env* env) {
  Assert(metrics != NULL);
  Assert(vectors != NULL);
  Assert(env != NULL);

  Insist(env, ( ! env->all2all ) ? "Unimplemented." : 0);

  /*---Denominator---*/

  Float_t* vector_sums = malloc(metrics->num_vector_local*sizeof(Float_t));

  compute_vector_sums(vectors, vector_sums, env);
 
  /*---Numerator---*/
 
  int i = 0;
  int j = 0;
  int k = 0;
  for (i = 0; i < metrics->num_vector_local; ++i) {
    for (j = i+1; j < metrics->num_vector_local; ++j) {
      for (k = j+1; k < metrics->num_vector_local; ++k)   {
        Float_t sum = 0;
        int field = 0;
        for (field = 0; field < vectors->num_field; ++field) {
          const Float_t value1 = Vectors_float_get(vectors, field, i, env);
          const Float_t value2 = Vectors_float_get(vectors, field, j, env);
          const Float_t value3 = Vectors_float_get(vectors, field, k, env);
          Float_t min12 = value1 < value2 ? value1 : value2;
          sum += min12;
          sum += value1 < value3 ? value1 : value3;
          sum += value2 < value3 ? value2 : value3;
          sum -= min12 < value3 ? min12 : value3;
        } /*---for field---*/
        Metrics_float_set_3(metrics, i, j, k, sum, env);
      } /*---for k---*/
    } /*---for j---*/
  } /*---for i---*/
  
  /*---Combine---*/

  for ( i = 0; i < metrics->num_vector_local; ++i ) {
    for ( j = i+1; j < metrics->num_vector_local; ++j ) {
      for ( k = j+1; k < metrics->num_vector_local; ++k) {
        const Float_t numerator = Metrics_float_get_3(metrics, i, j, k, env);
        const Float_t denominator = vector_sums[i] + vector_sums[j] + vector_sums[k];
        Metrics_float_set_3(metrics, i, j, k, 3 * numerator / (2 * denominator), env);
      } /*---for k---*/
    } /*---for j---*/
  } /*---for i---*/

  free( vector_sums );  

}

/*===========================================================================*/

void compute_metrics_czekanowski_3way_gpu(Metrics* metrics,
                                          Vectors* vectors,
                                          Env* env) {
  Assert(metrics != NULL);
  Assert(vectors != NULL);
  Assert(env != NULL);

  Insist(env, ( ! env->all2all ) ? "Unimplemented." : 0);

  /*---Denominator---*/

  Float_t* vector_sums = malloc(metrics->num_vector_local*sizeof(Float_t));

  compute_vector_sums(vectors, vector_sums, env);  

  /*---Numerator---*/

  int i = 0;
  int j = 0;
  int k = 0;
   
  const int numvec = metrics->num_vector_local;
  const int numfield = vectors->num_field;
 
  /*---Initialize MAGMA library---*/
  magma_minproduct_init();
    
  /*---Allocate magma CPU memory for vectors and for result---*/
  Float_t*  h_vectors = NULL; //Data matrix
  Float_t*  h_matM = NULL; //matrix matrix min product of X^T*X
  Float_t*  h_matV = NULL; //for fixed index j, column Vi = elementwise mins of Xj and Xi
  Float_t*  h_matB = NULL; //matrix matrix min product of V^T*X
#ifdef FP_PRECISION_DOUBLE
  magma_minproduct_dmalloc_pinned( &h_vectors, numvec*numfield);
  magma_minproduct_dmalloc_pinned( &h_matM,    numvec*numvec);
  magma_minproduct_dmalloc_pinned( &h_matV,    numfield*numvec);
  magma_minproduct_dmalloc_pinned( &h_matB,    numvec*numvec);
#endif
#ifdef FP_PRECISION_SINGLE
  magma_minproduct_smalloc_pinned( &h_vectors, numvec*numfield);
  magma_minproduct_smalloc_pinned( &h_matM,    numvec*numvec);
  magma_minproduct_smalloc_pinned( &h_matV,    numfield*numvec);
  magma_minproduct_smalloc_pinned( &h_matB,    numvec*numvec);
#endif
  
  /*---Allocate GPU mirrors for CPU arrays---*/
  Float_t* d_vectors = NULL;
  Float_t* d_matM = NULL;
  Float_t* d_matV = NULL;
  Float_t* d_matB = NULL;
#ifdef FP_PRECISION_DOUBLE
  magma_minproduct_dmalloc( &d_vectors, numvec*numfield);
  magma_minproduct_dmalloc( &d_matM,    numvec*numvec);
  magma_minproduct_dmalloc( &d_matV,    numfield*numvec);
  magma_minproduct_dmalloc( &d_matB,    numvec*numvec);
#endif
#ifdef FP_PRECISION_SINGLE
  magma_minproduct_smalloc( &d_vectors, numvec*numfield);
  magma_minproduct_smalloc( &d_matM,    numvec*numvec);
  magma_minproduct_smalloc( &d_matV,    numfield*numvec);
  magma_minproduct_smalloc( &d_matB,    numvec*numvec);
#endif

  /*---Copy in vectors---*/
  for (i = 0; i < numvec; ++i) {
    int k = 0;
    for (k = 0; k < numfield; ++k) {
      h_vectors[k+numfield*i] = Vectors_float_get( vectors, k, i, env );
    }
  }
  
  /*---Send matrix vectors to GPU---*/
#ifdef FP_PRECISION_DOUBLE  
  //magma_minproduct_dsetmatrix(numvec, numvec, h_matM, numvec,
  //                 d_matM, numvec);
  magma_minproduct_dsetmatrix(numfield, numvec, h_vectors, numfield,
                   d_vectors, numfield);
#endif
#ifdef FP_PRECISION_SINGLE
  //magma_minproduct_dsetmatrix(numvec, numvec, h_matM, numvec,
  //                 d_matM, numvec);
  magma_minproduct_dsetmatrix(numfield, numvec, h_vectors, numfield,
                   d_vectors, numfield);
#endif
    
  /*---Perform pseudo matrix-matrix product matM=vectors^T*vectors---*/
#ifdef FP_PRECISION_DOUBLE
  magma_minproductblas_dgemm
#endif
#ifdef FP_PRECISION_SINGLE
  magma_minproductblas_sgemm
#endif
   (Magma_minproductTrans,
    Magma_minproductNoTrans,
    numvec,
    numvec,
    numfield,
    1.0,
    d_vectors,
    numfield,
    d_vectors,
    numfield,
    0.0,
    d_matM,
    numvec); 
  
  /*---Copy matM from GPU---*/

#ifdef FP_PRECISION_DOUBLE
  magma_minproduct_dgetmatrix(numvec, numvec, d_matM, numvec, h_matM, numvec);
#endif
#ifdef FP_PRECISION_SINGLE
  magma_minproduct_sgetmatrix(numvec, numvec, d_matM, numvec, h_matM, numvec);
#endif

  /*---Process all combinations starting with j, i, k---*/
  for (j=1; j < numvec-1; ++j) {
    /*---Populate first j-1 columns of matV---*/
    for (i=0; i < numvec; ++i) {
      //Compare columns x_i and x_j element-wise
      for (k = 0; k < numfield; ++k) {
        //row = k, column = i
        h_matV[k+i*numfield] = h_vectors[k+numfield*i] < h_vectors[k+numfield*j] ? 
                             h_vectors[k+numfield*i] : h_vectors[k+numfield*j];
      }
    }

    /*---Send matrix matV to GPU---*/
//magma_minproduct_Xsetmatrix(M, N, host_matrix, leading_dim_of_host, device_matrix, leading_dim_of_dev);
#ifdef FP_PRECISION_DOUBLE
    magma_minproduct_dsetmatrix(numfield, numvec, h_matV, numfield, d_matV, numfield);
#endif
#ifdef FP_PRECISION_SINGLE
    magma_minproduct_ssetmatrix(numfield, numvec, h_matV, numfield, d_matV, numfield);
#endif
 
    /*---Perform matrix-matrix product matB = matV^T*vectors---*/
#ifdef FP_PRECISION_DOUBLE
    magma_minproductblas_dgemm
#endif
#ifdef FP_PRECISION_SINGLE
    magma_minproductblas_sgemm
#endif
    (Magma_minproductTrans, 
     Magma_minproductNoTrans,
     numvec, 
     numvec, 
     numfield, 
     1.0,
     d_matV, 
     numfield,  
     d_vectors, 
     numfield, 
     0.0,
     d_matB,
      numvec);
    
    /*---Copy matB from GPU---*/
#ifdef FP_PRECISION_DOUBLE
    magma_minproduct_dgetmatrix(numvec, numvec, d_matB, numvec, h_matB, numvec);
#endif
#ifdef FP_PRECISION_SINGLE
    magma_minproduct_sgetmatrix(numvec, numvec, d_matB, numvec, h_matB, numvec);
#endif
    /*---Combine results---*/
    for (i=0; i<j; ++i) {
      Float_t min_ij = h_matM[j+numvec*i];
      for (k = j+1; k < numvec; ++k) {
        Float_t min_ik = h_matM[k+numvec*i];
        Float_t min_jk = h_matM[k+numvec*j];
        //sum of mins vectors i, j, and k is matB(i,k)
        Float_t min_ijk = h_matB[k+numvec*i];
        const Float_t numerator = min_ij + min_ik + min_jk - min_ijk;
        const Float_t denominator = vector_sums[i] + vector_sums[j] + vector_sums[k];
        Metrics_float_set_3(metrics, i, j, k, 3 * numerator / (2 * denominator), env);
      }/*---for k---*/
    }/*---for i---*/
  }/*---for j---*/
    
  /*---Free memory and finalize---*/
  magma_minproduct_free(d_matM);
  magma_minproduct_free(d_vectors);
  magma_minproduct_free(d_matV);
  magma_minproduct_free(d_matB);
  magma_minproduct_free_pinned(h_matM);
  magma_minproduct_free_pinned(h_vectors);
  magma_minproduct_free_pinned(h_matB);
  magma_minproduct_free_pinned(h_matV);
    
  magma_minproduct_finalize();
    
  free(vector_sums);
}

