/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
       
       @author Azzam Haidar
       @author Stan Tomov
       @author Raffaele Solca
       
       @generated from zbulge_back_m.cpp normal z -> d, Fri Jan 30 19:00:18 2015

 */
#include "common_magma_tally2.h"
#include "magma_tally2_bulge.h"
#include "magma_tally2_dbulge.h"

#ifdef MAGMA_tally2_SETAFFINITY
#include "affinity.h"
#endif


#define PRECISION_d

static void *magma_tally2_dapplyQ_m_parallel_section(void *arg);

static void magma_tally2_dtile_bulge_applyQ(
    magma_tally2_int_t core_id, magma_tally2_side_t side, magma_tally2_int_t n_loc, magma_tally2_int_t n, magma_tally2_int_t nb, magma_tally2_int_t Vblksiz,
    double *E, magma_tally2_int_t lde,
    double *V, magma_tally2_int_t ldv,
    double *TAU,
    double *T, magma_tally2_int_t ldt);

//////////////////////////////////////////////////////////////////////////////////////////////////////////

class magma_tally2_dapplyQ_m_data {

public:

    magma_tally2_dapplyQ_m_data(magma_tally2_int_t ngpu_, magma_tally2_int_t threads_num_, magma_tally2_int_t n_, magma_tally2_int_t ne_, magma_tally2_int_t n_gpu_,
                         magma_tally2_int_t nb_, magma_tally2_int_t Vblksiz_, double *E_, magma_tally2_int_t lde_,
                         double *V_, magma_tally2_int_t ldv_, double *TAU_,
                         double *T_, magma_tally2_int_t ldt_)
    :
    ngpu(ngpu_),
    threads_num(threads_num_),
    n(n_),
    ne(ne_),
    n_gpu(n_gpu_),
    nb(nb_),
    Vblksiz(Vblksiz_),
    E(E_),
    lde(lde_),
    V(V_),
    ldv(ldv_),
    TAU(TAU_),
    T(T_),
    ldt(ldt_)
    {
        magma_tally2_int_t count = threads_num;

        if (threads_num > 1)
            --count;

        pthread_barrier_init(&barrier, NULL, count);
    }

    ~magma_tally2_dapplyQ_m_data()
    {
        pthread_barrier_destroy(&barrier);
    }
    const magma_tally2_int_t ngpu;
    const magma_tally2_int_t threads_num;
    const magma_tally2_int_t n;
    const magma_tally2_int_t ne;
    const magma_tally2_int_t n_gpu;
    const magma_tally2_int_t nb;
    const magma_tally2_int_t Vblksiz;
    double* const E;
    const magma_tally2_int_t lde;
    double* const V;
    const magma_tally2_int_t ldv;
    double* const TAU;
    double* const T;
    const magma_tally2_int_t ldt;
    pthread_barrier_t barrier;

private:

    magma_tally2_dapplyQ_m_data(magma_tally2_dapplyQ_m_data& data); // disable copy
};

class magma_tally2_dapplyQ_m_id_data {

public:

    magma_tally2_dapplyQ_m_id_data()
    : id(-1), data(NULL)
    {}

    magma_tally2_dapplyQ_m_id_data(magma_tally2_int_t id_, magma_tally2_dapplyQ_m_data* data_)
    : id(id_), data(data_)
    {}

    magma_tally2_int_t id;
    magma_tally2_dapplyQ_m_data* data;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////

extern "C" magma_tally2_int_t
magma_tally2_dbulge_back_m(
    magma_tally2_int_t ngpu,
    magma_tally2_uplo_t uplo,
    magma_tally2_int_t n, magma_tally2_int_t nb,
    magma_tally2_int_t ne, magma_tally2_int_t Vblksiz,
    double *Z, magma_tally2_int_t ldz,
    double *V, magma_tally2_int_t ldv,
    double *TAU,
    double *T, magma_tally2_int_t ldt,
    magma_tally2_int_t* info)
{
    magma_tally2_int_t threads = magma_tally2_get_parallel_numthreads();
    magma_tally2_int_t mklth   = magma_tally2_get_lapack_numthreads();
    magma_tally2_set_lapack_numthreads(1);

    real_Double_t timeaplQ2=0.0;

    double f= 1.;
    magma_tally2_int_t n_gpu = ne;

//#if defined(PRECISION_s) || defined(PRECISION_d)
//    double gpu_cpu_perf = 32; //gpu over cpu performance
//#else
//    double gpu_cpu_perf = 32;  // gpu over cpu performance
//#endif

    double perf_temp= .85;
    double perf_temp2= perf_temp;
    for (magma_tally2_int_t itmp=1; itmp < ngpu; ++itmp)
        perf_temp2 *= perf_temp;
    magma_tally2_int_t gpu_cpu_perf = magma_tally2_get_dbulge_gcperf();
    if (threads > 1) {
        f = 1. / (1. + (double)(threads-1)/ ((double)gpu_cpu_perf*(1.-perf_temp2)/(1.-perf_temp)));
        n_gpu = (magma_tally2_int_t)(f*ne);
    }






    /****************************************************
     *  apply V2 from left to the eigenvectors Z. dZ = (I-V2*T2*V2')*Z
     * **************************************************/

    timeaplQ2 = magma_tally2_wtime();

    /*============================
     *  use GPU+CPU's
     *==========================*/
//n_gpu = ne;
    if (n_gpu < ne) {
        // define the size of Q to be done on CPU's and the size on GPU's
        // note that GPU use Q(1:N_GPU) and CPU use Q(N_GPU+1:N)
        #ifdef ENABLE_DEBUG
        printf("---> calling GPU + CPU(if N_CPU > 0) to apply V2 to Z with NE %d     N_GPU %d   N_CPU %d\n",ne, n_gpu, ne-n_gpu);
        #endif
        magma_tally2_dapplyQ_m_data data_applyQ(ngpu, threads, n, ne, n_gpu, nb, Vblksiz, Z, ldz, V, ldv, TAU, T, ldt);

        magma_tally2_dapplyQ_m_id_data* arg;
        magma_tally2_malloc_cpu((void**) &arg, threads*sizeof(magma_tally2_dapplyQ_m_id_data));

        pthread_t* thread_id;
        magma_tally2_malloc_cpu((void**) &thread_id, threads*sizeof(pthread_t));

        pthread_attr_t thread_attr;

        // ===============================
        // relaunch thread to apply Q
        // ===============================
        // Set one thread per core
        pthread_attr_init(&thread_attr);
        pthread_attr_setscope(&thread_attr, PTHREAD_SCOPE_SYSTEM);
        pthread_setconcurrency(threads);

        // Launch threads
        for (magma_tally2_int_t thread = 1; thread < threads; thread++) {
            arg[thread] = magma_tally2_dapplyQ_m_id_data(thread, &data_applyQ);
            pthread_create(&thread_id[thread], &thread_attr, magma_tally2_dapplyQ_m_parallel_section, &arg[thread]);
        }
        arg[0] = magma_tally2_dapplyQ_m_id_data(0, &data_applyQ);
        magma_tally2_dapplyQ_m_parallel_section(&arg[0]);

        // Wait for completion
        for (magma_tally2_int_t thread = 1; thread < threads; thread++) {
            void *exitcodep;
            pthread_join(thread_id[thread], &exitcodep);
        }

        magma_tally2_free_cpu(thread_id);
        magma_tally2_free_cpu(arg);

        /*============================
         *  use only GPU
         *==========================*/
    } else {
        magma_tally2_dbulge_applyQ_v2_m(ngpu, Magma_tally2Left, ne, n, nb, Vblksiz, Z, ldz, V, ldv, T, ldt, info);
        magma_tally2_device_sync();
    }

    timeaplQ2 = magma_tally2_wtime()-timeaplQ2;

    magma_tally2_set_lapack_numthreads(mklth);
    return MAGMA_tally2_SUCCESS;
}

//##################################################################################################
static void *magma_tally2_dapplyQ_m_parallel_section(void *arg)
{
    magma_tally2_int_t my_core_id     = ((magma_tally2_dapplyQ_m_id_data*)arg) -> id;
    magma_tally2_dapplyQ_m_data* data = ((magma_tally2_dapplyQ_m_id_data*)arg) -> data;

    magma_tally2_int_t ngpu          = data -> ngpu;
    magma_tally2_int_t allcores_num   = data -> threads_num;
    magma_tally2_int_t n              = data -> n;
    magma_tally2_int_t ne             = data -> ne;
    magma_tally2_int_t n_gpu          = data -> n_gpu;
    magma_tally2_int_t nb             = data -> nb;
    magma_tally2_int_t Vblksiz        = data -> Vblksiz;
    double *E         = data -> E;
    magma_tally2_int_t lde            = data -> lde;
    double *V         = data -> V;
    magma_tally2_int_t ldv            = data -> ldv;
    double *TAU       = data -> TAU;
    double *T         = data -> T;
    magma_tally2_int_t ldt            = data -> ldt;
    pthread_barrier_t* barrier = &(data -> barrier);

    magma_tally2_int_t info;

    #ifdef ENABLE_TIMER
    real_Double_t timeQcpu=0.0, timeQgpu=0.0;
    #endif

    magma_tally2_int_t n_cpu = ne - n_gpu;

    // with MKL and when using omp_set_num_threads instead of mkl_set_num_threads
    // it need that all threads setting it to 1.
    magma_tally2_set_lapack_numthreads(1);

#ifdef MAGMA_tally2_SETAFFINITY
    //#define PRINTAFFINITY
#ifdef PRINTAFFINITY
    affinity_set print_set;
    print_set.print_affinity(my_core_id, "starting affinity");
#endif
    affinity_set original_set;
    affinity_set new_set(my_core_id);
    int check  = 0;
    int check2 = 0;
    // bind threads
    check = original_set.get_affinity();
    if (check == 0) {
        check2 = new_set.set_affinity();
        if (check2 != 0)
            printf("Error in sched_setaffinity (single cpu)\n");
    }
    else {
        printf("Error in sched_getaffinity\n");
    }
#ifdef PRINTAFFINITY
    print_set.print_affinity(my_core_id, "set affinity");
#endif
#endif

    if (my_core_id == 0) {
        //=============================================
        //   on GPU on thread 0:
        //    - apply V2*Z(:,1:N_GPU)
        //=============================================
        #ifdef ENABLE_TIMER
        timeQgpu = magma_tally2_wtime();
        #endif

        magma_tally2_dbulge_applyQ_v2_m(ngpu, Magma_tally2Left, n_gpu, n, nb, Vblksiz, E, lde, V, ldv, T, ldt, &info);
        magma_tally2_device_sync();

        #ifdef ENABLE_TIMER
        timeQgpu = magma_tally2_wtime()-timeQgpu;
        printf("  Finish Q2_GPU GGG timing= %f\n", timeQgpu);
        #endif
    } else {
        //=============================================
        //   on CPU on threads 1:allcores_num-1:
        //    - apply V2*Z(:,N_GPU+1:NE)
        //=============================================
        #ifdef ENABLE_TIMER
        if (my_core_id == 1)
            timeQcpu = magma_tally2_wtime();
        #endif

        magma_tally2_int_t n_loc = magma_tally2_ceildiv(n_cpu, allcores_num-1);
        double* E_loc = E + (n_gpu+ n_loc * (my_core_id-1))*lde;
        n_loc = min(n_loc,n_cpu - n_loc * (my_core_id-1));

        magma_tally2_dtile_bulge_applyQ(my_core_id, Magma_tally2Left, n_loc, n, nb, Vblksiz, E_loc, lde, V, ldv, TAU, T, ldt);
        pthread_barrier_wait(barrier);

        #ifdef ENABLE_TIMER
        if (my_core_id == 1) {
            timeQcpu = magma_tally2_wtime()-timeQcpu;
            printf("  Finish Q2_CPU CCC timing= %f\n", timeQcpu);
        }
        #endif
    } // END if my_core_id

#ifdef MAGMA_tally2_SETAFFINITY
    // unbind threads
    if (check == 0) {
        check2 = original_set.set_affinity();
        if (check2 != 0)
            printf("Error in sched_setaffinity (restore cpu list)\n");
    }
#ifdef PRINTAFFINITY
    print_set.print_affinity(my_core_id, "restored_affinity");
#endif
#endif

    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

#define E(m,n)   &(E[(m) + lde*(n)])
#define V(m)     &(V[(m)])
#define TAU(m)   &(TAU[(m)])
#define T(m)     &(T[(m)])
static void magma_tally2_dtile_bulge_applyQ(
    magma_tally2_int_t core_id, magma_tally2_side_t side, magma_tally2_int_t n_loc, magma_tally2_int_t n, magma_tally2_int_t nb, magma_tally2_int_t Vblksiz,
    double *E, magma_tally2_int_t lde,
    double *V, magma_tally2_int_t ldv,
    double *TAU,
    double *T, magma_tally2_int_t ldt)
    //, magma_tally2_int_t* info)
{
    //%===========================
    //%   local variables
    //%===========================
    magma_tally2_int_t firstcolj;
    magma_tally2_int_t bg, rownbm;
    magma_tally2_int_t st,ed,fst,vlen,vnb,colj;
    magma_tally2_int_t vpos,tpos;
    magma_tally2_int_t cur_blksiz,avai_blksiz, ncolinvolvd;
    magma_tally2_int_t nbgr, colst, coled;

    if (n <= 0)
        return;
    if (n_loc <= 0)
        return;

    //info = 0;
    magma_tally2_int_t INFO=0;

    magma_tally2_int_t nbGblk  = magma_tally2_ceildiv(n-1, Vblksiz);

    /*
     * version v1: for each chunck it apply all the V's then move to
     *                    the other chunck. the locality here inside each
     *                    chunck meaning that thread t apply V_k then move
     *                    to V_k+1 which overlap with V_k meaning that the
     *                    E_k+1 overlap with E_k. so here is the
     *                    locality however thread t had to read V_k+1 and
     *                    T_k+1 at each apply. note that all thread if they
     *                    run at same speed they might reading the same V_k
     *                    and T_k at the same time.
     * */

    magma_tally2_int_t nb_loc = 128; //$$$$$$$$

    magma_tally2_int_t     lwork = 2*nb_loc*max(Vblksiz,64);
    double *work, *work2;

    magma_tally2_dmalloc_cpu(&work, lwork);
    magma_tally2_dmalloc_cpu(&work2, lwork);

    magma_tally2_int_t nbchunk =  magma_tally2_ceildiv(n_loc, nb_loc);

    /* SIDE LEFT  meaning apply E = Q*E = (q_1*q_2*.....*q_n) * E ==> so traverse Vs in reverse order (forward) from q_n to q_1
     *            each q_i consist of applying V to a block of row E(row_i,:) and applies are overlapped meaning
     *            that q_i+1 overlap a portion of the E(row_i, :).
     *            IN parallel E is splitten in vertical block over the threads  */
    /* SIDE RIGHT meaning apply E = E*Q = E * (q_1*q_2*.....*q_n) ==> so tarverse Vs in normal  order (forward) from q_1 to q_n
     *            each q_i consist of applying V to a block of col E(:, col_i,:) and the applies are overlapped meaning
     *            that q_i+1 overlap a portion of the E(:, col_i).
     *            IN parallel E is splitten in horizontal block over the threads  */
    #ifdef ENABLE_DEBUG
    if ((core_id == 0) || (core_id == 1))
        printf("  APPLY Q2_cpu dbulge_back_m   N %d  N_loc %d  nbchunk %d  NB %d  Vblksiz %d  SIDE %c \n", n, n_loc, nbchunk, nb, Vblksiz, side);
    #endif
   
    for (magma_tally2_int_t i = 0; i < nbchunk; i++) {
        magma_tally2_int_t ib_loc = min(nb_loc, (n_loc - i*nb_loc));

        if (side == Magma_tally2Left) {
            for (bg = nbGblk; bg > 0; bg--) {
                firstcolj = (bg-1)*Vblksiz + 1;
                rownbm    = magma_tally2_ceildiv((n-(firstcolj+1)),nb);
                if (bg == nbGblk) rownbm    = magma_tally2_ceildiv((n-(firstcolj)),nb);  // last blk has size=1 used for real to handle A(N,N-1)
                for (magma_tally2_int_t j = rownbm; j > 0; j--) {
                    vlen = 0;
                    vnb  = 0;
                    colj      = (bg-1)*Vblksiz; // for k=0; I compute the fst and then can remove it from the loop
                    fst       = (rownbm -j)*nb+colj +1;
                    for (magma_tally2_int_t k=0; k < Vblksiz; k++) {
                        colj     = (bg-1)*Vblksiz + k;
                        st       = (rownbm -j)*nb+colj +1;
                        ed       = min(st+nb-1,n-1);
                        if (st > ed)
                            break;
                        if ((st == ed) && (colj != n-2))
                            break;
                        vlen=ed-fst+1;
                        vnb=k+1;
                    }
                    colst     = (bg-1)*Vblksiz;
                    magma_tally2_bulge_findVTpos(n, nb, Vblksiz, colst, fst, ldv, ldt, &vpos, &tpos);

                    if ((vlen > 0) && (vnb > 0)) {
                        lapackf77_dlarfb( "L", "N", "F", "C", &vlen, &ib_loc, &vnb, V(vpos), &ldv, T(tpos), &ldt, E(fst,i*nb_loc), &lde, work, &ib_loc);
                    }
                    if (INFO != 0)
                        printf("ERROR DORMQR INFO %d \n", (int) INFO);
                }
            }
        } else if (side == Magma_tally2Right) {
            rownbm    = magma_tally2_ceildiv((n-1),nb);
            for (magma_tally2_int_t k = 1; k <= rownbm; k++) {
                ncolinvolvd = min(n-1, k*nb);
                avai_blksiz=min(Vblksiz,ncolinvolvd);
                nbgr = magma_tally2_ceildiv(ncolinvolvd,avai_blksiz);
                for (magma_tally2_int_t j = 1; j <= nbgr; j++) {
                    vlen = 0;
                    vnb  = 0;
                    cur_blksiz = min(ncolinvolvd-(j-1)*avai_blksiz, avai_blksiz);
                    colst = (j-1)*avai_blksiz;
                    coled = colst + cur_blksiz -1;
                    fst   = (rownbm -k)*nb+colst +1;
                    for (colj=colst; colj <= coled; colj++) {
                        st       = (rownbm -k)*nb+colj +1;
                        ed       = min(st+nb-1,n-1);
                        if (st > ed)
                            break;
                        if ((st == ed) && (colj != n-2))
                            break;
                        vlen=ed-fst+1;
                        vnb=vnb+1;
                    }
                    magma_tally2_bulge_findVTpos(n, nb, Vblksiz, colst, fst, ldv, ldt, &vpos, &tpos);
                    if ((vlen > 0) && (vnb > 0)) {
                        lapackf77_dlarfb( "R", "N", "F", "C", &ib_loc, &vlen, &vnb, V(vpos), &ldv, T(tpos), &ldt, E(i*nb_loc,fst), &lde, work, &ib_loc);
                    }
                }
            }
        } else {
            printf("ERROR SIDE %d \n",side);
        }
    } // END loop over the chunks

    magma_tally2_free_cpu(work);
    magma_tally2_free_cpu(work2);
}
#undef E
#undef V
#undef TAU
#undef T
////////////////////////////////////////////////////////////////////////////////////////////////////
