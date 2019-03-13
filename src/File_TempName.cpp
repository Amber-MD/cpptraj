#include <cstdio>
#include <list>
#include "File_TempName.h"
#include "StringRoutines.h"
#include "CpptrajStdio.h"

namespace File {

typedef std::list<FileName> NameList;

static NameList TempFileNames_ = NameList();

// This seems a resonable limit for now
static const unsigned int maxtmpidx_ = 1628634;

static std::string TempPrefix_ = "cpptrajtmp";

/** Generate a temporary file name of format:
  *   TempPrefix_<#>
  * where <#> is based on the current number of temporary file names
  * that have been requested.
  */
FileName GenTempName() {
  // Could also set this to 0, but setting it to size is a better guarantee
  // that the name will be free.
  unsigned int tmpidx = TempFileNames_.size();
  FileName temp( TempPrefix_ + integerToString(tmpidx) );
  while (tmpidx < maxtmpidx_ && Exists(temp)) {
    tmpidx++;
    temp = FileName( TempPrefix_ + integerToString(tmpidx) );
  }
  if (tmpidx >= maxtmpidx_) {
    mprinterr("Internal Error: Too many temporary files. Remove files named '%s*'\n",
              TempPrefix_.c_str());
    return FileName();
  }
  return temp;
}

/** Free the given temporary file name. */
void FreeTempName(FileName const& temp) {
  TempFileNames_.remove( temp );
  if (Exists(temp)) remove(temp.full());
}
 
} /* END namespace File */
