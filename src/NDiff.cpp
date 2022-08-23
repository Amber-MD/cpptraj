#include "NDiff.h"
#include "BufferedLine.h"
#include "CpptrajStdio.h"

/** Intended to be a faster drop-in replacement for Nelson H. F. Beebe's 
  * ndiff.awk script.
  */
int NDiff(std::string const& fname1, std::string const& fname2, std::string const& tolarg,
          double tolIn)
{
  BufferedLine file1, file2;
  if (file1.OpenFileRead( fname1 )) {
    mprinterr("Error: ndiff: Could not open '%s'\n", fname1.c_str());
    return 1;
  }
  if (file2.OpenFileRead( fname2 )) {
    mprinterr("Error: ndiff: Could not open '%s'\n", fname2.c_str());
    return 1;
  }

  return 0;
}
