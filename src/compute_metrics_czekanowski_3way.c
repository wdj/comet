/****************************************************************/
/* Czekanowski Similarity Metric                                */
/*                                                              */
/* Code Author: Doug Hyatt                                      */
/* Created: June, 2015                                          */
/****************************************************************/

#include <sys/time.h>

#include "magma.h"
#include "magma_lapack.h"
#include <string.h>

#include "czek.h"

/*===========================================================================*/
long long binomial_coeff(int n, int k) {

  int i = 0;
  long long result;
  long long num = 1;
  long long den = 1;  

  if (n==k || k==0){
    result = 1;
  }
  else if (n < k){
    result = 0;
  }
  else {
    for (i=1; i < k+1; i++) {
      num *= (n + 1 - i);
      den *= i;
    }
    result = (long long) (num/den +.05);
  }
  return result;
}

/*===========================================================================*/

/* Function to map lexicographical combination to index */

long long index_of_comb(int n, int k, int *comb) {
  
  long long idx = 0;  
  int i = 0;
  int j = 0;
  int t[k+1];
  
  t[0]=0;
  for (i=0; i<k; ++i) {
    t[i+1]=comb[i]+1;
  }

  for (i = 1; i <= k; ++i) {
    if ((t[i-1]+1) <= (t[i]-1)) {
      for (j=t[i-1]+1; j <= (t[i]-1); ++j){
        idx += binomial_coeff(n-j, k-i);
      }
    }
  }
  return idx;
 
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

  long long counter = 0;
  long long n_choose_3 = binomial_coeff(numvec,3);
  czek_vals = malloc(n_choose_3*sizeof(Float_t));

  for (i = 0; i < metrics->num_vector_local; ++i) {
    for (j = i+1; j < metrics->num_vector_local; ++j) {
      for (k = j+1; k < metrics->num_vector_local; ++k)   {
        Float_t sum = 0;
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
        Metrics_Float_set_3(metrics, i, j, k, sum, env);
      } /*---for k---*/
    } /*---for j---*/
  } /*---for i---*/
  
  /*---Combine---*/

  for ( i = 0; i < metrics->num_vector_local; ++i ) {
    for ( j = i+1; j < metrics->num_vector_local; ++j ) {
      for ( k = j+1; k < metrics->num_vector_local; ++k) {
        const Float_t numerator = Metrics_Float_get_3(metrics, i, j, k, env);
        const Float_t denominator = vector_sums[i] + vector_sums[j] + vector_sums[k];
        Metrics_Float_set_3(metrics, i, j, k, 3 * numerator / (2 * denominator), env);
      } /*---for k---*/
    } /*---for j---*/
  } /*---for i---*/

  free( vector_sums );  

}

/*===========================================================================*/
void compute_metrics_czekanowski_3way_gpu(Metrics* metrics,
                                          Vectors* vectors,
                                          Env* env) {
    int i = 0;
    int j = 0;
    int k = 0;
    int k_index_offset = 0;
    int czek_comb[3];
    long long n_choose_3 = 0;
    long long czek_index = 0;
    Float_t* __restrict__ czek_vals = 0;
    Float_t* __restrict__ col_sums = 0;
    Float_t*   matX = 0; //Data matrix
    Float_t* d_matX = 0;
    Float_t*   matM = 0; //matrix matrix min product of X^T*X
    Float_t* d_matM = 0;
    Float_t*   matV = 0; //for fixed index j, column Vi = elementwise mins of Xj and Xi
    Float_t* d_matV = 0;
    Float_t*   matB = 0; //matrix matrix min product of V^T*X
    Float_t* d_matB = 0;
    Float_t min_ij = 0.0;
    Float_t min_ik = 0.0;
    Float_t min_jk = 0.0;
    Float_t min_ijk = 0.0;
    Float_t czek_value = 0.0;
    double time1 = 0.0;
    double time2 = 0.0;
    
    n_choose_3 = binomial_coeff(numvec,3);
    czek_vals = malloc(n_choose_3*sizeof(Float_t));
    col_sums = malloc(numvec*sizeof(Float_t));
    
    /* Initialize MAGMA library */
    magma_init();
    
    /* Allocate MAGMA CPU memory for vectors and for matM result */
    magma_dmalloc_pinned(&matX,numvec*numfield);
    magma_dmalloc_pinned(&matM, numvec*numvec);
    /* Copy in vectors to matrix*/
    for (j = 0; j < numvec; ++j) {
        for (i = 0; i < numfield; ++i) {
            matX[i+numfield*j] = vectors[j].data[i];
        }
    }
    /* Allocate GPU mirrors for CPU arrays */
    magma_dmalloc(&d_matX, numvec*numfield);
    magma_dmalloc(&d_matM, numvec*numvec);
    /* Initialize result to zero (apparently is required) */
    for (j = 0; j < numvec; ++j) {
        for (i = 0; i < numvec; ++i) {
            matM[i+numvec*j] = 0;
        }
    }
    
    /* Start timer */
    time1 = get_time();
    
    /* Compute individual vector sums */
    for (i = 0; i < numvec; ++i) {
        col_sums[i] = vector_sum(numfield, vectors[i].data);
    }
    
    /* Send matrices to GPU */
    magma_dsetmatrix(numvec, numvec, matM, numvec,
                     d_matM, numvec);
    magma_dsetmatrix(numfield, numvec, matX, numfield,
                     d_matX, numfield);
    
    /* Perform pseudo matrix-matrix product with data matrix*/
    magmablas_dgemm_tesla(MagmaTrans, MagmaNoTrans, numvec, 
                          numvec, numfield, 1.0, d_matX,
                          numfield, d_matX, numfield, 
                          0.0, d_matM, numvec);
    
    /* Copy result from GPU */
    magma_dgetmatrix(numvec, numvec, d_matM, numvec,
                     matM, numvec);

    /* Allocate MAGMA memory */
    //magma_dmalloc_pinned(&matV,         numfield*(numvec-1));
    magma_dmalloc_pinned(&matV,         numfield*numvec);
    //magma_dmalloc_pinned(&matB, (numvec-1)*numvec);
    magma_dmalloc_pinned(&matB, numvec*numvec);
    magma_dmalloc(&d_matV,          numfield*numvec);
    magma_dmalloc(&d_matB,   numvec*numvec);
    
    /* Initialize to zero */
    for (j=0; j<numvec; ++j) {
      for (i=0; i<numvec; ++i) {
        matB[i+numvec*j] = 0.0;
      }
    }
    for (j=0; j<numvec; ++j) {
      for (i=0; i<numfield; ++i) {
        matV[i+j*numfield] = 0.0; 
      }
    } 
    
    for (j=1; j < numvec-1; ++j) {
      /* Populate first j-1 columns of matrix A */
      for (i=0; i < numvec; ++i) {
        //Compare columns x_i and x_j element-wise
        for (k = 0; k < numfield; ++k) {
          //row = k, column = i
          matV[k+i*numfield] = min_op(matX[k+numfield*i], matX[k+numfield*j]);
         }
       }

        /* Send matrix to GPU */
        //magma_dsetmatrix(M, N, host_matrix, leading_dim_of_host, device_matrix, leading_dim_of_dev);
        magma_dsetmatrix(numfield, numvec, matV, numfield, d_matV, numfield);
        /* Perform matrix-matrix product */
        magmablas_dgemm_tesla(MagmaTrans, MagmaNoTrans,
                              numvec, numvec, numfield, 1.0,
                              d_matV, numfield, d_matX, numfield, 0.0,
                              d_matB, numvec);
        /* Copy result from GPU */
        magma_dgetmatrix(numvec, numvec, d_matB, numvec, matB, numvec);
        
        /* Compute 3-way Czek metric */
        //czek_comb[1]=j;
        for (i=0; i<j; ++i) {
            //czek_comb[0]=i;
            min_ij = matM[j+numvec*i];
            for (k = j+1; k < numvec; ++k) {
                //czek_comb[2]=k;
                //czek_index = index_of_comb(numvec,3,czek_comb);
                min_ik = matM[k+numvec*i];
                min_jk = matM[k+numvec*j];
                //comparison of vectors x_i, x_j, and x_k is matB(i,k)
                min_ijk = matB[k+numvec*i];
                /*
                if (matB[k+numvec*i]-matB[i+numvec*k]) {
                  printf("i%i j%i k%i, min_ijk = %f\n",i,j,k,min_ijk);
                }
                */
                czek_vals[czek_index++] = (((Float_t)3)*(min_ij + min_ik + min_jk - min_ijk)) /
                                       (((Float_t)2)*(col_sums[i] + col_sums[j] + col_sums[k]));
            }
        }
    }
    
    time2 = get_time();
///*
    Float_t* czek_vals_reorder = 0;
    czek_vals_reorder = malloc(n_choose_3*sizeof(Float_t));
    long long counter = 0;
    for (j=0; j<numvec-1; ++j) {
      czek_comb[1]=j;
      for (i=0; i<j; ++i) {
        czek_comb[0]=i; 
        for (k=j+1; k<numvec; ++k) { 
          czek_comb[2]=k;
          czek_index = index_of_comb(numvec,3,czek_comb);
          czek_vals_reorder[czek_index] = czek_vals[counter++];   
        }
      }
    }
       
     
  //  */

    /* Print output */ 
    czek_index=0;
    for (i = 0; i < numvec-2; ++i) {
        for (j = i+1; j < numvec-1; ++j) {
            for (k = j+1; k < numvec; ++k) {
                czek_value = czek_vals[czek_index++];
#ifndef NO_PRINT
                fprintf(fp, "%s\t%s\t%s\t%.4f\n", vectors[i].id, vectors[j].id, vectors[k].id, czek_value);
#endif
            }
        }
    }
    
    /* Print computation time and checksum */
    printf("alt2 numvec %i numfield %i "
           "time: %.6f "
           "checksum %.15e\n",
           numvec, numfield, time2-time1, checksum(czek_vals_reorder, numvec));
    
    
    /* Free memory */
    magma_free(d_matM);
    magma_free(d_matX);
    magma_free(d_matV);
    magma_free(d_matB);
    magma_free_pinned(matM);
    magma_free_pinned(matX);
    magma_free_pinned(matB);
    magma_free_pinned(matV);
    
    magma_finalize();
    
    free(czek_vals);
    free(col_sums);
    free(czek_vals_reorder);
    
}


/*===========================================================================*/

/* Read in the vectors (ids and data) and dynamically allocate */
/* memory as needed. */
struct _vector *read_vectors(FILE *fp, int *numvec, int *numfield) {
  int i = 0;
  int c = 0;               /* Character we read in */
  int field_ctr = 0;       /* Count the number of tab-delimited fields */
  int vec_ctr = 0;         /* Count the number of vectors (lines) */
  int ch_ctr = 0;          /* Counter for individual chars */
  char buf[MAXFIELDSIZE] = "";    /* Buffer for reading fields */
  char *conv_ptr = NULL;   /* Pointer to check conversion */
  struct _vector *vectors = NULL; /* Pointer to the vectors */

  int numvec_orig = 0;
  int numfield_orig = 0;

  /* Read in the first line to determine the number of fields */
  do {
    c = getc(fp);
    if (c == '\t') field_ctr++;
  } while(c != EOF && c != '\n');
  numfield_orig = field_ctr;
  *numfield = numfield_orig * NCOPIES_F;

  /* Now handle various error situations */
  if (feof(fp) != 0) {
    fprintf(stderr, "Read error: Saw end of file before expected.\n");
    return NULL;
  }
  if (ferror(fp) != 0 || c != '\n') {
    fprintf(stderr, "Read error: Error reading the file.\n");
    return NULL;
  }

  /* Allocate initial memory for the vectors and ids */
  vectors = (struct _vector *)malloc(VECCHUNK * sizeof(struct _vector) * NCOPIES_V);
  if (vectors == NULL) {
    fprintf(stderr, "Failed to allocate memory for the vectors.\n");
    return NULL;
  }
  for (i = 0; i < VECCHUNK; i++)
    memset(&vectors[i], 0, sizeof(struct _vector));

  /* Read in the remaining lines of the file and add data to the vectors */
  while(c != EOF) {

    /* Allocate memory for a new vector */
    vectors[vec_ctr].data = (Float_t *)malloc((*numfield)*sizeof(Float_t));
    if (vectors[vec_ctr].data == NULL) {
      fprintf(stderr, "Failed to allocate memory for the vectors.\n");
      free_vectors(vectors, vec_ctr);
      return NULL;
    }

    /* Reset field/character counters */
    field_ctr = 0;
    ch_ctr = 0;

    /* Read characters until we hit a tab or newline, then process */
    /* the buffer, either adding it to ID (field 1) or vector data. */
    do {
      c = getc(fp);
      if (c == '\t' || c == '\n') {
        if (field_ctr == 0) vectors[vec_ctr].id[ch_ctr] = '\0';
        else if (field_ctr <= numfield_orig) {
          buf[ch_ctr] = '\0';
          vectors[vec_ctr].data[field_ctr-1] = (Float_t)strtod(buf, &conv_ptr);
          if (conv_ptr != buf+ch_ctr) {
            fprintf(stderr, "Error converting '%s' to double, line %d.\n",
                    buf, vec_ctr+2);
            free_vectors(vectors, vec_ctr);
            return NULL;
          }
        }
        if (c == '\t') field_ctr++; 
        ch_ctr = 0;
      }
      else if (field_ctr == 0)
        vectors[vec_ctr].id[ch_ctr++] = c;
      else buf[ch_ctr++] = c;
    } while(c != EOF && c != '\n');
    if (c == EOF) break;

    /* If we see wrong number of fields... */
    if (field_ctr != numfield_orig) {
      fprintf(stderr, "Did not see correct number of fields.\n");
      fprintf(stderr, "Expected %d, saw %d, line %d.\n", numfield_orig, field_ctr, vec_ctr+2);
      free_vectors(vectors, vec_ctr);
      return NULL;
    }

    /* Increment vector counter.  If hit our limit, realloc for more. */
    vec_ctr++;
    if (vec_ctr%VECCHUNK == 0) {
      vectors = (struct _vector *)realloc(vectors,
                (vec_ctr+VECCHUNK) * sizeof(struct _vector) * NCOPIES_V);
      if (vectors == NULL) {
        free_vectors(vectors, vec_ctr);
        fprintf(stderr, "Failed to realloc memory for vectors.\n");
        return NULL;
      }
      for (i = vec_ctr; i < vec_ctr + VECCHUNK; i++)
        memset(&vectors[i], 0, sizeof(struct _vector));
    }
  }
  if (ferror(fp) != 0) {
    fprintf(stderr, "Read error: Error reading the file.\n");
    free_vectors(vectors, vec_ctr);
    return NULL;
  }

  numvec_orig = vec_ctr;
  *numvec = numvec_orig * NCOPIES_V;

  /* Replicate fields */
  for (i=0; i<numvec_orig; ++i) {
    int k = 0;
    for (k=0; k<numfield_orig; ++k) {
      int j = 0;
      for (j=1; j<NCOPIES_F; ++j) {
        vectors[i].data[k+j*numfield_orig] = vectors[i].data[k];
      }
    }
  }

  /* Replicate vectors */
  for (i=0; i<numvec_orig; ++i) {
    int j = 0;
    for (j=1; j<NCOPIES_V; ++j) {
      int k = 0;
      strcpy(vectors[i+j*numvec_orig].id, vectors[i].id);
      vectors[i+j*numvec_orig].data = (Float_t *)malloc((*numfield)*sizeof(Float_t));
      for (k=0; k<*numfield; ++k) {
        vectors[i+j*numvec_orig].data[k] = vectors[i].data[k];
      }
    }
  }

  return vectors;
}

/* Free routine for vectors */
void free_vectors(struct _vector *vectors, int numvec) {
  int i = 0;

  if (vectors != NULL) {
    for (i = 0; i < numvec; i++)
      if (vectors[i].data != NULL) free(vectors[i].data);
    free(vectors);
  } 
}
