/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c

*/
#include "common_magma_tally4.h"

#define A(m,n) (a+(n)*(*lda)+(m))

// task execution code
static void SCHED_zgemm(Quark* quark)
{
  magma_tally4_int_t UPPER;
  magma_tally4_int_t M;
  magma_tally4_int_t N;
  magma_tally4_int_t K;
  cuDoubleComplex *A1;
  magma_tally4_int_t LDA;
  cuDoubleComplex *A2;
  cuDoubleComplex *A3;

  cuDoubleComplex mone = MAGMA_tally4_Z_NEG_ONE;
  cuDoubleComplex one = MAGMA_tally4_Z_ONE;

  quark_unpack_args_8(quark, UPPER, M, N, K, A1, LDA, A2, A3);

  if (UPPER) {

    blasf77_zgemm(Magma_tally4ConjTransStr, Magma_tally4NoTransStr, 
      &M, &N, &K, &mone, A1, &LDA, A2, &LDA, &one, A3, &LDA);

  } else {

    blasf77_zgemm(Magma_tally4NoTransStr, Magma_tally4ConjTransStr, 
      &M, &N, &K, &mone, A1, &LDA, A2, &LDA, &one, A3, &LDA);

  }

}

// task execution code
void SCHED_zsyrk(Quark* quark)
{
  magma_tally4_int_t UPPER;
  magma_tally4_int_t N;
  cuDoubleComplex *A1;
  magma_tally4_int_t LDA;
  magma_tally4_int_t K;
  cuDoubleComplex *A2;

  double mone = -1.0;
  double one = 1.0;

  quark_unpack_args_6(quark, UPPER, N, K, A1, LDA, A2);

  if (UPPER) {

    blasf77_zherk(Magma_tally4UpperStr, Magma_tally4ConjTransStr, &N, &K, &mone, A1, &LDA, &one, 
      A2, &LDA);

  } else {

    blasf77_zherk(Magma_tally4LowerStr, Magma_tally4NoTransStr, &N, &K, &mone, A1, &LDA, &one, 
      A2, &LDA);

  }

}

// task execution code
void SCHED_zpotrf(Quark* quark)
{
  magma_tally4_int_t UPPER;
  magma_tally4_int_t N;
  cuDoubleComplex *A;
  magma_tally4_int_t LDA;

  magma_tally4_int_t *iinfo;

  magma_tally4_int_t info;

  quark_unpack_args_5(quark, UPPER, N, A, LDA, iinfo);

  if (UPPER) {

    lapackf77_zpotrf(Magma_tally4UpperStr, &N, A, &LDA, &info);

  } else {

    lapackf77_zpotrf(Magma_tally4LowerStr, &N, A, &LDA, &info);

  }

  if (info > 0) {
    iinfo[1] = iinfo[0] + info;
  }

}

// task execution code
void SCHED_ztrsm(Quark* quark)
{
  magma_tally4_int_t UPPER;
  magma_tally4_int_t M;
  magma_tally4_int_t N;
  cuDoubleComplex *A1;
  magma_tally4_int_t LDA;
  cuDoubleComplex *A2;

  cuDoubleComplex one = MAGMA_tally4_Z_ONE;

  quark_unpack_args_6(quark, UPPER, M, N, A1, LDA, A2);

  if (UPPER) {

    blasf77_ztrsm(Magma_tally4LeftStr, Magma_tally4UpperStr, Magma_tally4ConjTransStr, Magma_tally4NonUnitStr,
      &M, &N, &one, A1, &LDA, A2, &LDA);

  } else {

    blasf77_ztrsm(Magma_tally4RightStr, Magma_tally4LowerStr, Magma_tally4ConjTransStr, Magma_tally4NonUnitStr,
      &M, &N, &one, A1, &LDA, A2, &LDA);

  }

}

extern "C" magma_tally4_int_t 
magma_tally4_zpotrf_mc(magma_tally4_context *cntxt, char *uplo,
        magma_tally4_int_t *n,
        cuDoubleComplex *a, magma_tally4_int_t *lda,
        magma_tally4_int_t *info)
{
/*  -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

    Purpose   
    =======   
    ZPOTRF computes the Cholesky factorization of a Hermitian   
    positive definite matrix A.   

    The factorization has the form   
       A = U**T * U,  if UPLO = 'U', or   
       A = L  * L**T,  if UPLO = 'L',   
    where U is an upper triangular matrix and L is lower triangular.   

    This is the block version of the algorithm, calling Level 3 BLAS.   

    Arguments   
    =========   
    CNTXT   (input) MAGMA_tally4_CONTEXT
            CNTXT specifies the MAGMA_tally4 hardware context for this routine.   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) COMPLEX_16 array, dimension (LDA,N)   
            On entry, the Hermitian matrix A.  If UPLO = 'U', the leading   
            N-by-N upper triangular part of A contains the upper   
            triangular part of the matrix A, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading N-by-N lower triangular part of A contains the lower   
            triangular part of the matrix A, and the strictly upper   
            triangular part of A is not referenced.   

            On exit, if INFO = 0, the factor U or L from the Cholesky   
            factorization A = U**T*U or A = L*L**T.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, the leading minor of order i is not   
                  positive definite, and the factorization could not be   
                  completed.   
    =====================================================================   */

  if (cntxt->num_cores == 1 && cntxt->num_gpus == 1)
  {
    //magma_tally4_int_t result = magma_tally4_zpotrf(*uplo, *n, a, *lda, info);
    //return result;
  }

  // check arguments
  magma_tally4_int_t upper = (magma_tally4_int_t) lsame_(uplo, "U");                                          
  *info = 0;
  if (! upper && ! lsame_(uplo, "L")) {
    *info = -1;
  } else if (*n < 0) {
    *info = -2;
  } else if (*lda < max(1,*n)) {
    *info = -4;
  }
    if (*info != 0) {
        magma_tally4_xerbla( __func__, -(*info) );
        return MAGMA_tally4_ERR_ILLEGAL_VALUE;
    }

  Quark* quark = cntxt->quark;

  // get block size
  magma_tally4_int_t nb = (cntxt->nb ==-1)? magma_tally4_get_zpotrf_nb(*n): cntxt->nb;

  magma_tally4_int_t i,j,k;
  magma_tally4_int_t ii,jj,kk;
  magma_tally4_int_t temp,temp2,temp3;

  char label[10000];

  magma_tally4_int_t iinfo[2];
  iinfo[1] = 0;
  ii = -1;

  // traverse diagonal blocks
  for (i = 0; i < (*n); i += nb) {
    ii++;
    temp2 = min(nb,(*n)-i);
 
    // if not first block
    if (i > 0) {

      // first do large syrk, then split
      if (i < (*n)/2) {

        sprintf(label, "SYRK %d", ii);

        if (upper) {

          QUARK_Insert_Task(quark, SCHED_zsyrk, 0,
            sizeof(magma_tally4_int_t),             &upper,    VALUE,
            sizeof(magma_tally4_int_t),             &temp2,    VALUE,
            sizeof(magma_tally4_int_t),             &i,        VALUE,
            sizeof(cuDoubleComplex)*(*n)*(*n), A(0,i),    INPUT,
            sizeof(magma_tally4_int_t),             lda,       VALUE,
            sizeof(cuDoubleComplex)*(*n)*(*n), A(i,i),    INOUT,
            sizeof(cuDoubleComplex)*(*n)*(*n), A(i-nb,i), INPUT,
            strlen(label)+1,         label,     VALUE | TASKLABEL,
            6,                       "green",   VALUE | TASKCOLOR,
            0);

        } else {

          QUARK_Insert_Task(quark, SCHED_zsyrk, 0,
            sizeof(magma_tally4_int_t),             &upper,    VALUE,
            sizeof(magma_tally4_int_t),             &temp2,    VALUE,
            sizeof(magma_tally4_int_t),             &i,        VALUE,
            sizeof(cuDoubleComplex)*(*n)*(*n), A(i,0),    INPUT,
            sizeof(magma_tally4_int_t),             lda,       VALUE,
            sizeof(cuDoubleComplex)*(*n)*(*n), A(i,i),    INOUT,
            sizeof(cuDoubleComplex)*(*n)*(*n), A(i,i-nb), INPUT,
            strlen(label)+1,         label,     VALUE | TASKLABEL,
            6,                       "green",   VALUE | TASKCOLOR,
            0);
        }   

      } else {

        jj = -1;

        // split syrk into tiles
        for (j = 0; j < i; j += nb) {
          jj++;

          sprintf(label, "SYRK %d %d", ii, jj);

          if (upper) {

            QUARK_Insert_Task(quark, SCHED_zsyrk, 0,
              sizeof(magma_tally4_int_t),             &upper,    VALUE,
              sizeof(magma_tally4_int_t),             &temp2,    VALUE,
              sizeof(magma_tally4_int_t),             &nb,       VALUE,
              sizeof(cuDoubleComplex)*(*n)*(*n), A(j,i),    INPUT,
              sizeof(magma_tally4_int_t),             lda,       VALUE,
              sizeof(cuDoubleComplex)*(*n)*(*n), A(i,i),    INOUT,
              strlen(label)+1,         label,     VALUE | TASKLABEL,
              6,                       "green",   VALUE | TASKCOLOR,
              0);

          } else {

            QUARK_Insert_Task(quark, SCHED_zsyrk, 0,
              sizeof(magma_tally4_int_t),             &upper,    VALUE,
              sizeof(magma_tally4_int_t),             &temp2,    VALUE,
              sizeof(magma_tally4_int_t),             &nb,       VALUE,
              sizeof(cuDoubleComplex)*(*n)*(*n), A(i,j),    INPUT,
              sizeof(magma_tally4_int_t),             lda,       VALUE,
              sizeof(cuDoubleComplex)*(*n)*(*n), A(i,i),    INOUT,
              strlen(label)+1,         label,     VALUE | TASKLABEL,
              6,                       "green",   VALUE | TASKCOLOR,
              0);
          }
        }
      }

      // if not last block
      if (i < ((*n)-nb)) {

        jj = -1;

        // split gemm into tiles
        for (j = i+nb; j < (*n); j += nb){
      jj++;
          kk = -1;

          for (k = 0; k < i; k += nb) {
            kk++;
            temp = min(nb,(*n)-j);

            sprintf(label, "GEMM %d %d %d", ii, jj, kk);

            if (upper) {
              QUARK_Insert_Task(quark, SCHED_zgemm, 0,
                sizeof(magma_tally4_int_t),             &upper,    VALUE,
                sizeof(magma_tally4_int_t),             &nb,       VALUE,
                sizeof(magma_tally4_int_t),             &temp,     VALUE,
                sizeof(magma_tally4_int_t),             &nb,       VALUE,
                sizeof(cuDoubleComplex)*(*n)*(*n), A(k,i), INPUT,
                sizeof(magma_tally4_int_t),             lda,       VALUE,
                sizeof(cuDoubleComplex)*(*n)*(*n), A(k,j),    INPUT,
                sizeof(cuDoubleComplex)*(*n)*(*n), A(i,j), INOUT,
                strlen(label)+1,         label,     VALUE | TASKLABEL,
                5,                       "blue",    VALUE | TASKCOLOR,
                0);

            } else {

              QUARK_Insert_Task(quark, SCHED_zgemm, 0,
                sizeof(magma_tally4_int_t),             &upper,    VALUE,
                sizeof(magma_tally4_int_t),             &temp,     VALUE,
                sizeof(magma_tally4_int_t),             &nb,       VALUE,
                sizeof(magma_tally4_int_t),             &nb,       VALUE,
                sizeof(cuDoubleComplex)*(*n)*(*n), A(j,k), INPUT,
                sizeof(magma_tally4_int_t),             lda,       VALUE,
                sizeof(cuDoubleComplex)*(*n)*(*n), A(i,k),    INPUT,
                sizeof(cuDoubleComplex)*(*n)*(*n), A(j,i), INOUT,
                strlen(label)+1,         label,     VALUE | TASKLABEL,
                5,                       "blue",    VALUE | TASKCOLOR,
                0);
            }
          }
        }
      }  
    }

    iinfo[0] = i;

    sprintf(label, "POTRF %d", ii);

    QUARK_Insert_Task(quark, SCHED_zpotrf, 0,
      sizeof(magma_tally4_int_t),             &upper,    VALUE,
      sizeof(magma_tally4_int_t),             &temp2,    VALUE,
      sizeof(cuDoubleComplex)*(*n)*(*n), A(i,i),    INOUT,
      sizeof(magma_tally4_int_t),             lda,       VALUE,
      sizeof(magma_tally4_int_t),             iinfo,     OUTPUT,
      strlen(label)+1,         label,     VALUE | TASKLABEL,
      5,                       "cyan",    VALUE | TASKCOLOR,
      0);

    // if not last block
    if (i < ((*n)-nb)) {

      // split trsm into tiles
      for (j = i + nb; j < (*n); j += nb) {

        temp = min(nb,(*n)-j);

        sprintf(label, "TRSM %d", ii);

        if (upper) {

          QUARK_Insert_Task(quark, SCHED_ztrsm, 0,
            sizeof(magma_tally4_int_t),             &upper,    VALUE,
            sizeof(magma_tally4_int_t),             &nb,       VALUE,
            sizeof(magma_tally4_int_t),             &temp,     VALUE,
            sizeof(cuDoubleComplex)*(*n)*(*n), A(i,i),    INPUT,
            sizeof(magma_tally4_int_t),             lda,       VALUE,
            sizeof(cuDoubleComplex)*(*n)*(*n), A(i,j),    INOUT,
            strlen(label)+1,         label,     VALUE | TASKLABEL,
            4,                       "red",     VALUE | TASKCOLOR,
            0);

        } else {

          QUARK_Insert_Task(quark, SCHED_ztrsm, 0,
            sizeof(magma_tally4_int_t),             &upper,    VALUE,
            sizeof(magma_tally4_int_t),             &temp,     VALUE,
            sizeof(magma_tally4_int_t),             &nb,       VALUE,
            sizeof(cuDoubleComplex)*(*n)*(*n), A(i,i),    INPUT,
            sizeof(magma_tally4_int_t),             lda,       VALUE,
            sizeof(cuDoubleComplex)*(*n)*(*n), A(j,i),    INOUT,
            strlen(label)+1,         label,     VALUE | TASKLABEL,
            4,                       "red",     VALUE | TASKCOLOR,
            0);
        }
      } 
    }
  }
  
  QUARK_Barrier(quark);
}

#undef A

