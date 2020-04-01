//-----------------------------------------------------------------------------
/*!
 * \file   input_output.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  I/O functions used by driver, header.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _comet_vectors_io_hh_
#define _comet_vectors_io_hh_

#include "env.hh"
#include "vectors.hh"

//=============================================================================

namespace comet {

//=============================================================================
/// \brief Perform input/output operations for Vectors class.

class VectorsIO {
public:

  void static read(GMVectors& vectors, const char* path, CEnv& env);
  void static write(GMVectors& vectors, const char* path, CEnv& env);

};

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _comet_vectors_io_hh_

//-----------------------------------------------------------------------------
