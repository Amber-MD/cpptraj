#ifndef FILE_TEMPNAME_H
#define FILE_TEMPNAME_H
#include "FileName.h"
namespace File {
  /// \return A temporary file name.
  FileName GenTempName();
  /// Indicate temporary file name no longer needed.
  void FreeTempName(FileName const&);
} /* END namespace File */
#endif
