//-----------------------------------------------------------------------------
/*!
 * \file   assertions.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Macros and code for assertions.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _gm_assertions_hh_
#define _gm_assertions_hh_

//=============================================================================
// Macros.

//-----------------------------------------------------------------------------
/// \brief Insist macro - assert a condition even for release builds.
///
///        This should be used only for non-performance-sensitive
///        code locations -- e.g., not in a deep loop nest.

#define GMInsist(condition) \
  (void)((condition) || (CoMet::assert_(#condition, __FILE__, __LINE__), 0))

//-----------------------------------------------------------------------------
/// \brief Insist macro specifically for a user-caused error condition.

#define GMInsistInterface(env, condition) \
  (void)((condition) || \
         (CoMet::insist_interface(env, #condition, __FILE__, __LINE__), 0))

//-----------------------------------------------------------------------------
/// \brief Assertion macro (for debug builds only).

#ifndef NDEBUG
#define GM_ASSERTIONS_ON
#define GMAssert(condition) GMInsist(condition)
#else
#define GMAssert(condition)
#endif

//-----------------------------------------------------------------------------
/// \brief Static (i.e., compile time) assertion macro.

// TODO: replace with C++11 equivalent.

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

//=============================================================================
// Declarations.

namespace CoMet {

//-----------------------------------------------------------------------------
/// \brief Function to support the GMAssert, GMInsist macros.

//        Use trailing underscore to avoid collision with C assert macro.

void assert_(const char* condition_string, const char* file, int line);

//-----------------------------------------------------------------------------
/// \brief Function to support the GMInsistInterface macro.

void insist_interface(void const * const env,
                         const char* condition_string,
                         const char* file,
                         int line);

//=============================================================================

} // namespace CoMet

//-----------------------------------------------------------------------------

#endif // _gm_assertions_hh_

//-----------------------------------------------------------------------------
