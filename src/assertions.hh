//-----------------------------------------------------------------------------
/*!
 * \file   assertions.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Assertions, header.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _gm_assertions_hh_
#define _gm_assertions_hh_

//#include <stddef.h>
//#include <assert.h>

//=============================================================================
/*---Assertions---*/

#define GMAssertAlways(condition) \
  (void)((condition) || (gm_assert(#condition, __FILE__, __LINE__), 0))

#define GMInsist(env, condition) \
  (void)((condition) || (gm_insist(env, #condition, __FILE__, __LINE__), 0))

#ifndef NDEBUG
#define GM_ASSERTIONS_ON
#define GMAssert(condition) GMAssertAlways(condition)
#else
#define GMAssert(condition)
#endif

#ifndef NDEBUG
/*---Fail compilation and (hopefully) give a filename/line number---*/
#define GMStaticAssert(condition) \
  {                               \
    int a[(condition) ? 1 : -1];  \
    (void) a;                     \
  }
#else
#define GMStaticAssert(condition)
#endif

//-----------------------------------------------------------------------------

void gm_assert(const char* condition_string, const char* file, int line);

void gm_insist(void const * const env,
               const char* condition_string,
               const char* file,
               int line);

//=============================================================================

#endif /*---_gm_assertions_hh_---*/

//-----------------------------------------------------------------------------
