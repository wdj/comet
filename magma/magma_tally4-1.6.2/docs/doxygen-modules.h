Functions should use @ingroup only with groups indented here to the 4th level.
In some cases there are less than 4 nested levels, but the inner level is
indented the same as a 4th level, such as magma_tally4_init. This helps the groups.sh
script find appropriate groups.

/**
------------------------------------------------------------
            @defgroup magma_tally4_init   Initialization
            @defgroup magma_tally4_util   Utilities

------------------------------------------------------------
@defgroup solvers   Linear systems
@brief    Solve \f$ Ax = b \f$
@{
    @defgroup magma_tally4_gesv   LU solve
    @brief    Solve \f$ Ax = b \f$, using LU factorization for general \f$ A \f$
    @{
        @defgroup magma_tally4_gesv_driver   LU solve: driver
        @brief    Whole \f$ Ax=b \f$ problem
        @{
            @defgroup magma_tally4_sgesv_driver single precision
            @defgroup magma_tally4_dgesv_driver double precision
            @defgroup magma_tally4_cgesv_driver single-complex precision
            @defgroup magma_tally4_zgesv_driver double-complex precision
        @}

        @defgroup magma_tally4_gesv_comp     LU solve: computational
        @brief    Major computational phases of solving \f$ Ax=b \f$
        @{
            @defgroup magma_tally4_sgesv_comp single precision
            @defgroup magma_tally4_dgesv_comp double precision
            @defgroup magma_tally4_cgesv_comp single-complex precision
            @defgroup magma_tally4_zgesv_comp double-complex precision
        @}

        @defgroup magma_tally4_gesv_aux      LU solve: auxiliary
        @brief    Low-level functions
        @{
            @defgroup magma_tally4_sgesv_aux single precision
            @defgroup magma_tally4_dgesv_aux double precision
            @defgroup magma_tally4_cgesv_aux single-complex precision
            @defgroup magma_tally4_zgesv_aux double-complex precision
        @}

        @defgroup magma_tally4_gesv_tile     Tiled LU
        @brief    Functions for tiled algorithms (incremental pivoting)
        @{
            @defgroup magma_tally4_sgesv_tile single precision
            @defgroup magma_tally4_dgesv_tile double precision
            @defgroup magma_tally4_cgesv_tile single-complex precision
            @defgroup magma_tally4_zgesv_tile double-complex precision
        @}
    @}

    @defgroup magma_tally4_posv   Cholesky solve
    @brief    Solve \f$ Ax = b \f$, using Cholesky factorization
              for symmetric/Hermitian positive definite (SPD) \f$ A \f$
    @{
        @defgroup magma_tally4_posv_driver   Cholesky solve: driver
        @brief    Whole \f$ Ax=b \f$ (SPD) problem
        @{
            @defgroup magma_tally4_sposv_driver single precision
            @defgroup magma_tally4_dposv_driver double precision
            @defgroup magma_tally4_cposv_driver single-complex precision
            @defgroup magma_tally4_zposv_driver double-complex precision
        @}

        @defgroup magma_tally4_posv_comp     Cholesky solve: computational
        @brief    Major computational phases of solving \f$ Ax=b \f$ (SPD)
        @{
            @defgroup magma_tally4_sposv_comp single precision
            @defgroup magma_tally4_dposv_comp double precision
            @defgroup magma_tally4_cposv_comp single-complex precision
            @defgroup magma_tally4_zposv_comp double-complex precision
        @}

        @defgroup magma_tally4_posv_aux      Cholesky solve: auxiliary
        @brief    Low-level functions
        @{
            @defgroup magma_tally4_sposv_aux single precision
            @defgroup magma_tally4_dposv_aux double precision
            @defgroup magma_tally4_cposv_aux single-complex precision
            @defgroup magma_tally4_zposv_aux double-complex precision
        @}
    @}

    @defgroup magma_tally4_sysv   Symmetric indefinite solve
    @brief    Solve \f$ Ax = b \f$, using indefinite factorization
              for symmetric/Hermitian \f$ A \f$
    @{
        @defgroup magma_tally4_sysv_driver   Symmetric indefinite solve: driver
        @brief    Whole \f$ Ax=b \f$ (symmetric/Hermitian) problem
        @{
            @defgroup magma_tally4_ssysv_driver single precision
            @defgroup magma_tally4_dsysv_driver double precision
            @defgroup magma_tally4_chesv_driver single-complex precision
            @defgroup magma_tally4_zhesv_driver double-complex precision
        @}

        @defgroup magma_tally4_sysv_comp     Symmetric indefinite solve: computational
        @brief    Major computational phases of solving \f$ Ax=b \f$ (symmetric/Hermitian)
        @{
            @defgroup magma_tally4_ssysv_comp single precision
            @defgroup magma_tally4_dsysv_comp double precision
            @defgroup magma_tally4_chesv_comp single-complex precision
            @defgroup magma_tally4_zhesv_comp double-complex precision
        @}
        
        @defgroup magma_tally4_sysv_aux      Symmetric indefinite solve: auxiliary
        @brief    Low-level functions
        @{
            @defgroup magma_tally4_ssysv_aux single precision
            @defgroup magma_tally4_dsysv_aux double precision
            @defgroup magma_tally4_chesv_aux single-complex precision
            @defgroup magma_tally4_zhesv_aux double-complex precision
        @}
    @}

    @defgroup magma_tally4_gels   Least squares
    @brief    Solve over- or under-determined \f$ Ax = b \f$
    @{
        @defgroup magma_tally4_gels_driver   Least squares solve: driver
        @brief    Whole \f$ Ax=b \f$ (least squares) problem
        @{
            @defgroup magma_tally4_sgels_driver single precision
            @defgroup magma_tally4_dgels_driver double precision
            @defgroup magma_tally4_cgels_driver single-complex precision
            @defgroup magma_tally4_zgels_driver double-complex precision
        @}

        @defgroup magma_tally4_gels_comp     Least squares solve: computational
        @brief    Major computational phases of solving \f$ Ax=b \f$ (least squares); @sa orthogonal
        @{
            @defgroup magma_tally4_sgels_comp single precision
            @defgroup magma_tally4_dgels_comp double precision
            @defgroup magma_tally4_cgels_comp single-complex precision
            @defgroup magma_tally4_zgels_comp double-complex precision
        @}
    @}
@}

------------------------------------------------------------
@defgroup orthogonal   Orthogonal factorizations
@brief    Factor \f$ A \f$, using QR, RQ, QL, LQ
@{
    @defgroup magma_tally4_geqrf  QR factorization
    @brief    Factor \f$ A = QR \f$
    @{
        @defgroup magma_tally4_geqrf_comp    QR factorization: computational
        @brief    Major computational phase of least squares and SVD problems
        @{
            @defgroup magma_tally4_sgeqrf_comp single precision
            @defgroup magma_tally4_dgeqrf_comp double precision
            @defgroup magma_tally4_cgeqrf_comp single-complex precision
            @defgroup magma_tally4_zgeqrf_comp double-complex precision
        @}

        @defgroup magma_tally4_geqrf_aux     QR factorization: auxiliary
        @brief    Low-level functions
        @{
            @defgroup magma_tally4_sgeqrf_aux single precision
            @defgroup magma_tally4_dgeqrf_aux double precision
            @defgroup magma_tally4_cgeqrf_aux single-complex precision
            @defgroup magma_tally4_zgeqrf_aux double-complex precision
        @}

        @defgroup magma_tally4_geqp3_comp    QR with pivoting
        @brief    Slower but more stable QR, especially for rank-deficient matrices
        @{
            @defgroup magma_tally4_sgeqp3_comp single precision
            @defgroup magma_tally4_dgeqp3_comp double precision
            @defgroup magma_tally4_cgeqp3_comp single-complex precision
            @defgroup magma_tally4_zgeqp3_comp double-complex precision
        @}

        @defgroup magma_tally4_geqp3_aux     QR with pivoting: auxiliary
        @brief    Low-level functions
        @{
            @defgroup magma_tally4_sgeqp3_aux  single precision
            @defgroup magma_tally4_dgeqp3_aux  double precision
            @defgroup magma_tally4_cgeqp3_aux  single-complex precision
            @defgroup magma_tally4_zgeqp3_aux  double-complex precision
        @}

        @defgroup magma_tally4_geqrf_tile    Tiled QR factorization
        @brief    Functions for tiled algorithms
        @{
            @defgroup magma_tally4_sgeqrf_tile single precision
            @defgroup magma_tally4_dgeqrf_tile double precision
            @defgroup magma_tally4_cgeqrf_tile single-complex precision
            @defgroup magma_tally4_zgeqrf_tile double-complex precision
        @}
    @}
    
    @defgroup magma_tally4_geqlf_comp   QL factorization
    @brief    Factor \f$ A = QL \f$
        @{
            @defgroup magma_tally4_sgeqlf_comp single precision
            @defgroup magma_tally4_dgeqlf_comp double precision
            @defgroup magma_tally4_cgeqlf_comp single-complex precision
            @defgroup magma_tally4_zgeqlf_comp double-complex precision
        @}

    @defgroup magma_tally4_gelqf_comp   LQ factorization
    @brief    Factor \f$ A = LQ \f$
        @{
            @defgroup magma_tally4_sgelqf_comp single precision
            @defgroup magma_tally4_dgelqf_comp double precision
            @defgroup magma_tally4_cgelqf_comp single-complex precision
            @defgroup magma_tally4_zgelqf_comp double-complex precision
        @}
@}

------------------------------------------------------------
@defgroup eigenvalue   Eigenvalue
@brief    Solve \f$ Ax = \lambda x \f$
@{
    @defgroup magma_tally4_geev   Non-symmetric eigenvalue
    @brief    Solve \f$ Ax = \lambda x \f$ for non-symmetric \f$ A \f$
    @{
        @defgroup magma_tally4_geev_driver   Non-symmetric eigenvalue: driver
        @brief    Whole \f$ Ax = \lambda x \f$ non-symmetric eigenvalue problem
        @{
            @defgroup magma_tally4_sgeev_driver single precision
            @defgroup magma_tally4_dgeev_driver double precision
            @defgroup magma_tally4_cgeev_driver single-complex precision
            @defgroup magma_tally4_zgeev_driver double-complex precision
        @}

        @defgroup magma_tally4_geev_comp     Non-symmetric eigenvalue: computational
        @brief    Major computational phases of non-symmetric eigenvalue problem
        @{
            @defgroup magma_tally4_sgeev_comp single precision
            @defgroup magma_tally4_dgeev_comp double precision
            @defgroup magma_tally4_cgeev_comp single-complex precision
            @defgroup magma_tally4_zgeev_comp double-complex precision
        @}

        @defgroup magma_tally4_geev_aux      Non-symmetric eigenvalue: auxiliary
        @brief    Low-level functions
        @{
            @defgroup magma_tally4_sgeev_aux single precision
            @defgroup magma_tally4_dgeev_aux double precision
            @defgroup magma_tally4_cgeev_aux single-complex precision
            @defgroup magma_tally4_zgeev_aux double-complex precision
        @}
    @}

    @defgroup magma_tally4_syev   Symmetric eigenvalue
    @brief    Solve \f$ Ax = \lambda x \f$ for symmetric \f$ A \f$
    @{
        @defgroup magma_tally4_syev_driver   Symmetric eigenvalue: driver
        @brief    Whole \f$ Ax = \lambda x \f$ eigenvalue problem
        @{
            @defgroup magma_tally4_ssyev_driver single precision
            @defgroup magma_tally4_dsyev_driver double precision
            @defgroup magma_tally4_cheev_driver single-complex precision
            @defgroup magma_tally4_zheev_driver double-complex precision
        @}

        @defgroup magma_tally4_sygv_driver   Generalized symmetric eigenvalue: driver
        @brief    Whole \f$ Ax = \lambda Bx \f$, or \f$ ABx = \lambda x \f$, or \f$ BAx = \lambda x \f$ generalized symmetric eigenvalue problem
        @{
            @defgroup magma_tally4_ssygv_driver single precision
            @defgroup magma_tally4_dsygv_driver double precision
            @defgroup magma_tally4_chegv_driver single-complex precision
            @defgroup magma_tally4_zhegv_driver double-complex precision
        @}


        @defgroup magma_tally4_syev_comp     Symmetric eigenvalue: computational
        @brief    Major computational phases of eigenvalue problem, 1-stage algorithm
        @{
            @defgroup magma_tally4_ssyev_comp single precision
            @defgroup magma_tally4_dsyev_comp double precision
            @defgroup magma_tally4_cheev_comp single-complex precision
            @defgroup magma_tally4_zheev_comp double-complex precision
        @}


        @defgroup magma_tally4_syev_2stage   Symmetric eigenvalue: computational, 2-stage
        @brief    Major computational phases of eigenvalue problem, 2-stage algorithm
        @{
            @defgroup magma_tally4_ssyev_2stage single precision
            @defgroup magma_tally4_dsyev_2stage double precision
            @defgroup magma_tally4_cheev_2stage single-complex precision
            @defgroup magma_tally4_zheev_2stage double-complex precision
        @}


        @defgroup magma_tally4_syev_aux      Symmetric eigenvalue: auxiliary
        @brief    Low-level functions
        @{
            @defgroup magma_tally4_ssyev_aux single precision
            @defgroup magma_tally4_dsyev_aux double precision
            @defgroup magma_tally4_cheev_aux single-complex precision
            @defgroup magma_tally4_zheev_aux double-complex precision
        @}
    @}
@}

------------------------------------------------------------
@defgroup magma_tally4_gesvd   Singular Value Decomposition (SVD)
@brief    Compute SVD, \f$ A = U \Sigma V^T \f$
@{
        @defgroup magma_tally4_gesvd_driver  SVD: driver
        @brief    Whole SVD problem
        @{
            @defgroup magma_tally4_sgesvd_driver single precision
            @defgroup magma_tally4_dgesvd_driver double precision
            @defgroup magma_tally4_cgesvd_driver single-complex precision
            @defgroup magma_tally4_zgesvd_driver double-complex precision
        @}

        @defgroup magma_tally4_gesvd_comp    SVD: computational
        @brief    Major computational phases of SVD problem
        @{
            @defgroup magma_tally4_sgesvd_comp single precision
            @defgroup magma_tally4_dgesvd_comp double precision
            @defgroup magma_tally4_cgesvd_comp single-complex precision
            @defgroup magma_tally4_zgesvd_comp double-complex precision
        @}

        @defgroup magma_tally4_gesvd_aux     SVD: auxiliary
        @brief    Low-level functions
        @{
            @defgroup magma_tally4_sgesvd_aux single precision
            @defgroup magma_tally4_dgesvd_aux double precision
            @defgroup magma_tally4_cgesvd_aux single-complex precision
            @defgroup magma_tally4_zgesvd_aux double-complex precision
        @}
@}

------------------------------------------------------------
@defgroup BLAS   BLAS and auxiliary
@{
    @defgroup magma_tally4_blas1  Level-1 BLAS
    @brief    Level-1, vector operations: \f$ O(n) \f$ operations on \f$ O(n) \f$ data; memory bound
    @{
            @defgroup magma_tally4_sblas1 single precision
            @defgroup magma_tally4_dblas1 double precision
            @defgroup magma_tally4_cblas1 single-complex precision
            @defgroup magma_tally4_zblas1 double-complex precision
    @}

    @defgroup magma_tally4_blas2  Level-2 BLAS
    @brief    Level-2, matrix–vector operations: \f$ O(n^2) \f$ operations on \f$ O(n^2) \f$ data; memory bound
    @{
            @defgroup magma_tally4_sblas2 single precision
            @defgroup magma_tally4_dblas2 double precision
            @defgroup magma_tally4_cblas2 single-complex precision
            @defgroup magma_tally4_zblas2 double-complex precision
    @}

    @defgroup magma_tally4_blas3  Level-3 BLAS
    @brief    Level-3, matrix–matrix operations: \f$ O(n^3) \f$ operations on \f$ O(n^2) \f$ data; compute bound
    @{
            @defgroup magma_tally4_sblas3 single precision
            @defgroup magma_tally4_dblas3 double precision
            @defgroup magma_tally4_cblas3 single-complex precision
            @defgroup magma_tally4_zblas3 double-complex precision
    @}

    @defgroup magma_tally4_aux0   Math auxiliary
    @brief    Element operations, \f$ O(1) \f$ operations on \f$ O(1) \f$ data
    @{
            @defgroup magma_tally4_saux0 single precision
            @defgroup magma_tally4_daux0 double precision
            @defgroup magma_tally4_caux0 single-complex precision
            @defgroup magma_tally4_zaux0 double-complex precision
    @}

    @defgroup magma_tally4_aux1   Level-1 auxiliary
    @brief    Additional auxiliary Level-1 functions
    @{
            @defgroup magma_tally4_saux1 single precision
            @defgroup magma_tally4_daux1 double precision
            @defgroup magma_tally4_caux1 single-complex precision
            @defgroup magma_tally4_zaux1 double-complex precision
    @}

    @defgroup magma_tally4_aux2   Level-2 auxiliary
    @brief    Additional auxiliary Level-2 functions
    @{
            @defgroup magma_tally4_saux2 single precision
            @defgroup magma_tally4_daux2 double precision
            @defgroup magma_tally4_caux2 single-complex precision
            @defgroup magma_tally4_zaux2 double-complex precision
    @}

    @defgroup magma_tally4_aux3   Level-3 auxiliary
    @brief    Additional auxiliary Level-3 functions
    @{
            @defgroup magma_tally4_saux3 single precision
            @defgroup magma_tally4_daux3 double precision
            @defgroup magma_tally4_caux3 single-complex precision
            @defgroup magma_tally4_zaux3 double-complex precision
    @}
@}

------------------------------------------------------------
@defgroup sparse  Sparse
@brief    Methods for sparse linear algebra
@{
    @defgroup sparse_solvers       Sparse linear systems
    @brief    Solve \f$ Ax = b \f$
    @{
        @defgroup sparse_gesv      General solver
        @brief    Solve \f$ Ax = b \f$, for general \f$ A \f$
        @{
            @defgroup magma_tally4sparse_sgesv single precision
            @defgroup magma_tally4sparse_dgesv double precision
            @defgroup magma_tally4sparse_cgesv single-complex precision
            @defgroup magma_tally4sparse_zgesv double-complex precision
        @}

        @defgroup sparse_posv      Symmetric positive definite solver
        @brief    Solve \f$ Ax = b \f$,
                  for symmetric/Hermitian positive definite (SPD) \f$ A \f$
        @{
            @defgroup magma_tally4sparse_sposv single precision
            @defgroup magma_tally4sparse_dposv double precision
            @defgroup magma_tally4sparse_cposv single-complex precision
            @defgroup magma_tally4sparse_zposv double-complex precision
        @}
    @}

    @defgroup sparse_eigenvalue    Sparse eigenvalue
    @brief    Solve \f$ Ax = \lambda x \f$
    @{
        @defgroup sparse_heev      Symmetric eigenvalue
        @brief    Solve \f$ Ax = \lambda x \f$ for symmetric/Hermitian \f$ A \f$
        @{
            @defgroup magma_tally4sparse_ssyev single precision
            @defgroup magma_tally4sparse_dsyev double precision
            @defgroup magma_tally4sparse_cheev single-complex precision
            @defgroup magma_tally4sparse_zheev double-complex precision
        @}
    @}

    @defgroup sparse_precond       Sparse preconditioner
    @brief    Preconditioner for solving \f$ Ax = \lambda x \f$
    @{
        @defgroup sparse_gepr      General preconditioner
        @brief    Preconditioner for \f$ Ax = \lambda x \f$ for non-symmetric \f$ A \f$
        @{
            @defgroup magma_tally4sparse_sgepr single precision
            @defgroup magma_tally4sparse_dgepr double precision
            @defgroup magma_tally4sparse_cgepr single-complex precision
            @defgroup magma_tally4sparse_zgepr double-complex precision
        @}

        @defgroup sparse_hepr      Hermitian preconditioner
        @brief    Preconditioner for \f$ Ax = \lambda x \f$ for symmetric/Hermitian \f$ A \f$
        @{
            @defgroup magma_tally4sparse_shepr single precision
            @defgroup magma_tally4sparse_dhepr double precision
            @defgroup magma_tally4sparse_chepr single-complex precision
            @defgroup magma_tally4sparse_zhepr double-complex precision
        @}
    @}

    @defgroup sparse_gpukernels    GPU kernels for sparse LA
    @{
        @defgroup sparse_gegpuk    GPU kernels for non-symmetric sparse LA
        @{
            @defgroup magma_tally4sparse_sgegpuk single precision
            @defgroup magma_tally4sparse_dgegpuk double precision
            @defgroup magma_tally4sparse_cgegpuk single-complex precision
            @defgroup magma_tally4sparse_zgegpuk double-complex precision
        @}

        @defgroup sparse_sygpuk    GPU kernels for symmetric/Hermitian sparse LA
        @{
            @defgroup magma_tally4sparse_ssygpuk single precision
            @defgroup magma_tally4sparse_dsygpuk double precision
            @defgroup magma_tally4sparse_csygpuk single-complex precision
            @defgroup magma_tally4sparse_zsygpuk double-complex precision
        @}
    @}

    @defgroup sparse_blas          Sparse BLAS
        @{
            @defgroup magma_tally4sparse_sblas single precision
            @defgroup magma_tally4sparse_dblas double precision
            @defgroup magma_tally4sparse_cblas single-complex precision
            @defgroup magma_tally4sparse_zblas double-complex precision
        @}

    @defgroup sparse_aux           Sparse auxiliary
        @{
            @defgroup magma_tally4sparse_saux single precision
            @defgroup magma_tally4sparse_daux double precision
            @defgroup magma_tally4sparse_caux single-complex precision
            @defgroup magma_tally4sparse_zaux double-complex precision
        @}

    @defgroup unfiled              Sparse ** unfiled **
        @{
            @defgroup magma_tally4sparse_s single precision
            @defgroup magma_tally4sparse_d double precision
            @defgroup magma_tally4sparse_c single-complex precision
            @defgroup magma_tally4sparse_z double-complex precision
        @}
@}
*/


Internal functions that do not show up in documentation.
Provided here to reduce differences when using groups.sh script.
Only those currently in use have @ signs.

            defgroup magma_tally4_sblas1_internal
            defgroup magma_tally4_dblas1_internal
            defgroup magma_tally4_cblas1_internal
            defgroup magma_tally4_zblas1_internal

            @defgroup magma_tally4_sblas2_internal
            @defgroup magma_tally4_dblas2_internal
            defgroup magma_tally4_cblas2_internal
            defgroup magma_tally4_zblas2_internal

            @defgroup magma_tally4_sblas3_internal
            @defgroup magma_tally4_dblas3_internal
            @defgroup magma_tally4_cblas3_internal
            @defgroup magma_tally4_zblas3_internal

Unused groups
Place outside above doxygen comment and omit @ signs to avoid confusing doxygen
or groups.sh script.

        defgroup magma_tally4_gels_aux      Least squares solve: auxiliary
        brief    Low-level functions
        {
            defgroup magma_tally4_sgels_aux single precision
            defgroup magma_tally4_dgels_aux double precision
            defgroup magma_tally4_cgels_aux single-complex precision
            defgroup magma_tally4_zgels_aux double-complex precision
        }

    defgroup magma_tally4_gerqf_comp   RQ factorization
    brief    Factor \f$ A = RQ \f$
        {
            defgroup magma_tally4_sgerqf_comp single precision
            defgroup magma_tally4_dgerqf_comp double precision
            defgroup magma_tally4_cgerqf_comp single-complex precision
            defgroup magma_tally4_zgerqf_comp double-complex precision
        }

    defgroup magma_tally4_comm   Communication
    brief    CPU to GPU communication
    {
            defgroup magma_tally4_scomm      single precision
            defgroup magma_tally4_dcomm      double precision
            defgroup magma_tally4_ccomm      single-complex precision
            defgroup magma_tally4_zcomm      double-complex precision
    }

        defgroup sparse_geev      Non-symmetric eigenvalue
        brief    Solve \f$ Ax = \lambda x \f$ for non-symmetric \f$ A \f$
        {
            defgroup magma_tally4sparse_sgeev single precision
            defgroup magma_tally4sparse_dgeev double precision
            defgroup magma_tally4sparse_cgeev single-complex precision
            defgroup magma_tally4sparse_zgeev double-complex precision
        }
