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

//-----------------------------------------------------------------------------

// Check for error from improper user settings

#define GMInsistInterface(env, condition) \
  (void)((condition) || \
         (gm_insist_interface(env, #condition, __FILE__, __LINE__), 0))

// Insist - assert always, incl. release mode

#define GMInsist(condition) \
  (void)((condition) || (gm_assert(#condition, __FILE__, __LINE__), 0))

// Assertion

#ifndef NDEBUG
#define GM_ASSERTIONS_ON
#define GMAssert(condition) GMInsist(condition)
#else
#define GMAssert(condition)
#endif

// Static assertion

#ifndef NDEBUG
// Fail compilation and (hopefully) give a filename/line number
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

void gm_insist_interface(void const * const env,
                         const char* condition_string,
                         const char* file,
                         int line);

//=============================================================================

#endif // _gm_assertions_hh_

//-----------------------------------------------------------------------------
