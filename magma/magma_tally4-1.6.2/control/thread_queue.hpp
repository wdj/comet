/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Mark Gates
*/

#ifndef MAGMA_tally4_THREAD_HPP
#define MAGMA_tally4_THREAD_HPP

#include <queue>

#include "common_magma_tally4.h"


// ---------------------------------------------
extern "C"
void* magma_tally4_thread_main( void* arg );


// ---------------------------------------------
class magma_tally4_task
{
public:
    magma_tally4_task() {}
    virtual ~magma_tally4_task() {}
    
    virtual void run() = 0;  // pure virtual function to execute task
};


// ---------------------------------------------
// Thread pool with multi-producer, multi-consumer queue.
//
// This is similar to python's queue class, but also implements worker threads
// and adds quit mechanism.
// sync is like python's join. Threads do not exit, so I find join to be a misleading name.
class magma_tally4_thread_queue
{
public:
    magma_tally4_thread_queue();
    ~magma_tally4_thread_queue();
    
    void launch( magma_tally4_int_t in_nthread );
    void push_task( magma_tally4_task* task );
    void sync();
    void quit();
    
protected:
    friend void* magma_tally4_thread_main( void* arg );
    magma_tally4_task* pop_task();
    void task_done();
    
    magma_tally4_int_t get_thread_index( pthread_t thread ) const;
    
private:
    std::queue< magma_tally4_task* > q;  ///<  queue of tasks
    bool            quit_flag;    ///<  quit() sets this to true; after this, pop returns NULL
    magma_tally4_int_t     ntask;        ///<  number of unfinished tasks (in queue or currently executing)
    pthread_mutex_t mutex;        ///<  mutex lock for queue, quit, ntask
    pthread_cond_t  cond;         ///<  condition variable for changes to queue and quit (see push, pop, quit)
    pthread_cond_t  cond_ntask;   ///<  condition variable for changes to ntask (see sync, task_done)
    pthread_t*      threads;      ///<  array of threads
    magma_tally4_int_t     nthread;      ///<  number of threads
};

#endif        //  #ifndef MAGMA_tally4_THREAD_HPP
