/**
 *
 * magma_tally4winthread.h
 *
 *  This file handles the mapping from pthreads calls to windows threads.
 *  MAGMA_tally4 is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.3.1
 * @author Piotr Luszczek
 * @date January 2015
 *
 * This file is originally from PLASMA project, where plasma has been
 * replaced by MAGMA_tally4.
 *
 **/
#ifndef MAGMA_tally4WINTHREAD_H
#define MAGMA_tally4WINTHREAD_H

#if (defined( _WIN32 ) || defined( _WIN64 )) && ! defined( __MINGW32__ )

#include <windows.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct pthread_s {
  HANDLE hThread;
  unsigned int uThId;
} pthread_t;

typedef HANDLE pthread_mutex_t;
typedef int pthread_mutexattr_t;
typedef int pthread_attr_t;
typedef int pthread_condattr_t;

typedef struct pthread_cond_s {
  HANDLE hSem;
  HANDLE hEvt;
  CRITICAL_SECTION cs;
  int waitCount; /* waiting thread counter */
} pthread_cond_t;

typedef int pthread_attr_t;

#define PTHREAD_MUTEX_INITIALIZER ((pthread_mutex_t) -1)

#define PTHREAD_SCOPE_SYSTEM 1

#define MAGMA_tally4_DLLPORT
#define MAGMA_tally4_CDECL __cdecl

MAGMA_tally4_DLLPORT pthread_t MAGMA_tally4_CDECL pthread_self(void);
MAGMA_tally4_DLLPORT int MAGMA_tally4_CDECL pthread_mutex_init(pthread_mutex_t *mutex, const pthread_mutexattr_t * attr);
MAGMA_tally4_DLLPORT int MAGMA_tally4_CDECL pthread_mutex_destroy(pthread_mutex_t *mutex);
MAGMA_tally4_DLLPORT int MAGMA_tally4_CDECL pthread_mutex_lock(pthread_mutex_t *mutex);
MAGMA_tally4_DLLPORT int MAGMA_tally4_CDECL pthread_mutex_trylock(pthread_mutex_t *mutex);
MAGMA_tally4_DLLPORT int MAGMA_tally4_CDECL pthread_mutex_unlock(pthread_mutex_t *mutex);
MAGMA_tally4_DLLPORT int MAGMA_tally4_CDECL pthread_attr_init(pthread_attr_t *attr);
MAGMA_tally4_DLLPORT int MAGMA_tally4_CDECL pthread_attr_destroy(pthread_attr_t *attr);
MAGMA_tally4_DLLPORT int MAGMA_tally4_CDECL pthread_attr_setscope(pthread_attr_t *attr, int scope);
MAGMA_tally4_DLLPORT int MAGMA_tally4_CDECL pthread_create(pthread_t *tid, const pthread_attr_t *attr, void *(*start) (void *), void *arg);
MAGMA_tally4_DLLPORT int MAGMA_tally4_CDECL pthread_cond_init(pthread_cond_t *cond, const pthread_condattr_t *attr);
MAGMA_tally4_DLLPORT int MAGMA_tally4_CDECL pthread_cond_destroy(pthread_cond_t *cond);
MAGMA_tally4_DLLPORT int MAGMA_tally4_CDECL pthread_cond_wait(pthread_cond_t *cond, pthread_mutex_t *mutex);
MAGMA_tally4_DLLPORT int MAGMA_tally4_CDECL pthread_cond_broadcast(pthread_cond_t *cond);
MAGMA_tally4_DLLPORT int MAGMA_tally4_CDECL pthread_join(pthread_t thread, void **value_ptr);
MAGMA_tally4_DLLPORT int MAGMA_tally4_CDECL pthread_equal(pthread_t thread1, pthread_t thread2);

MAGMA_tally4_DLLPORT int MAGMA_tally4_CDECL pthread_setconcurrency (int);

MAGMA_tally4_DLLPORT unsigned int MAGMA_tally4_CDECL pthread_self_id(void);

#ifdef __cplusplus
}
#endif

#endif /* (_WIN32 || _WIN64) && ! __MINGW32__ */

#endif /* MAGMA_tally4WINTHREAD_H */