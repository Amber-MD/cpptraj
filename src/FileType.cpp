#include "FileType.h"
#include <sys/stat.h> // stat

Cpptraj::File::Type Cpptraj::File::IdType(std::string const& filenameIn) {
  if (filenameIn.empty()) return T_UNKNOWN;

  struct stat frame_stat;
  if (stat(filenameIn.c_str(), &frame_stat) == -1) {
    return T_UNKNOWN;
  }

  if ( frame_stat.st_mode & S_IFDIR ) {
    return T_DIR;
  } else if ( frame_stat.st_mode & S_IFREG ) {
    return T_FILE;
  }
  return T_UNKNOWN;
}
